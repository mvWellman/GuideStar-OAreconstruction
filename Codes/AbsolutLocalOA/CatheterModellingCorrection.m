function [LHS,RHS] = CatheterModellingCorrection(M1,M2)
% CatheterModellingCorrection Computes the LHS and RHS correction matrices
% for accurate optic axis measurements.
%
%   This function is used for correcting the cumulative round-trip signals
%   from the inside (M1) and outside (M2) of a birefringent sheath, which
%   is assumed to be a Q-oriented linear retarder. The function derives the
%   LHS and RHS multipliers needed to obtain absolute optic axis
%   measurements from these signals. Input:
%       M1 - Cumulative round-trip signal from the inside of the
%       birefringent sheath. M2 - Cumulative round-trip signal from the
%       outside of the birefringent sheath.
%
%   Output:
%       LHS - The left-hand side correction matrix. RHS - The right-hand
%       side correction matrix.


% find Lsys
% Decompose rotation matrices R1 into their components
m1 = decomposeRot(M1);
% Define the function to find the optimal L1 linear retarder
fun = @(x) findLsys(x,M1);
% Initial guess for the optimization (mean of the decomposed rotation)
start = double(mean(m1,2)/2);

% Optimize to find the best L1 parameters
xopt = fminsearch(fun,start([1;2]));

% Compute the corrected rotation matrices using the optimal parameters
[~,m1C] = fun(xopt);

% Create rotation matrices for the correction
Lsys = makeRot3x3(cat(1,-xopt,0));
inverseRootSquareM1 = makeRot3x3(m1C);

% Unwrap the phase angle for the corrected rotation matrices
phi = unwrap(atan2(m1C(2,:),m1C(1,:)));

% negative rotation: Adjust the rotation matrix if the mean difference in phi is negative
if mean(diff(phi))<0
    Lsys = makeRot3x3(cat(1,(-xopt+pi*xopt/norm(xopt)),0));
    inverseRootSquareM1 = pagemtimes(Lsys,pagemtimes(M1,Lsys));
end


% just to make sure we have no wrap, unwrap this signal (in most instances
% this is unnecessary)
m1C = unwrapOA1D(inverseRootSquareM1);
inverseRootSquareM1 = makeRot3x3(-m1C/2);

% correct R2 with the same Linv
inverseRootSquareM2 = pagemtimes(Lsys,pagemtimes(M2,Lsys));

% Compute the sheath signal
isolatedSheath = pagemtimes(inverseRootSquareM1,pagemtimes(inverseRootSquareM2,inverseRootSquareM1));
sheathSignal = unwrapOA1D(isolatedSheath);

% also compute the 2pi wrapped version
isolatedSheath2 = sheathSignal - 2*pi*sheathSignal./sqrt(sum(sheathSignal.^2,1));

%pick the one with the smaller ret
if mean(sqrt(sum(sheathSignal.^2,1)))>mean(sqrt(sum(isolatedSheath2.^2,1)))
    sheathSignal = isolatedSheath2;
end

% Create a rotation matrix to align with the Q-axis
Vtheta = makeRot3x3(cat(1,zeros(2,size(M1,3)),-atan2(m1C(2,:),m1C(1,:))));

% Rotate the sheath signal with Vtheta
rloc = squeeze(pagemtimes(Vtheta,permute(sheathSignal,[1,3,2])));

% we want to find the theta offset that aligns the rotation vectors with
% the Q-axis. Considering SO3 matrices, this corresponds to matrices that
% are of the form [1,0,0;0,m22,m23;0,m32,m33]. Hence, we reduce the energy
% of the elements that we want to be zero:
fun = @(x)squeeze(sum(sum((makeRot(pagemtimes(makeRot3x3(cat(1,zeros(2,numel(x)),x)),permute(rloc,[1,3,4,2])))).^2.*[0;1;1;1;0;0;1;0;0],1),4));
angle1 = fminsearch(fun,0);
en1 = fun(angle1);

% this could be a local minimum; check optimim +pi/2
angle2 = fminsearch(fun,angle1+pi/2);
en2 = fun(angle2);%

if en2<en1
    meanAngle = angle2;
else
    meanAngle = angle1;
end

% Update the rotation matrix with the optimal angle
Vtheta = makeRot3x3(cat(1,zeros(2,size(M1,3)),-atan2(m1C(2,:),m1C(1,:))+meanAngle));
rloc2 = squeeze(pagemtimes(Vtheta,permute(sheathSignal,[1,3,2])));

% make sure this is along positive Q
if mean(rloc2(1,:))<0
    meanAngle = meanAngle + pi;
    Vtheta = makeRot3x3(cat(1,zeros(2,size(M1,3)),-atan2(m1C(2,:),m1C(1,:))+meanAngle));
end

%
RHS = pagemtimes(Lsys,pagemtimes(inverseRootSquareM1,pagemtimes(makeRot3x3(-sheathSignal/2),pagetranspose(Vtheta))));
LHS = pagetranspose(RHS.*[1,1,-1;1,1,-1;-1,-1,1]);

%% 

function [err,oacenter] = findLsys(params,Mcath)
% [err,oacenter] = findLsys(params,Mball) computes the error of the lsq trace
% minimization problem. oacenter is the optic axis of the catheter sheath
% signal after correcting for the linear element L1.

% Compensates Mcath with L1inv, which is a linear retarder
L1inv = makeRot3x3(-[params(1);params(2);0]);
Mcomp = pagemtimes(L1inv,pagemtimes(Mcath,L1inv));

% compute unwrapped ret and then compute the trace minimum
oacenter = decomposeRot(Mcomp); % unwrapping is necessary, as this is an artifact of the varying optic axis orientation; 
ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(oacenter,[],2).^2,1))>pi,2),2)),[3,1])>0;
delta = 2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
oacenter(ss) = oacenter(ss) - delta(ss);
ret = sqrt(sum(oacenter.^2,1));
if mean(ret)>pi % the modulo 2*pi solutions do not change the sense of orientation, and simply switch the sign of mret
    oacenter = oacenter-2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
    ret = sqrt(sum(oacenter.^2,1));
end
mret = atan2(mean(sin(real(ret))),mean(cos(real(ret))));
err = 1-mean(cos(ret-mret));
