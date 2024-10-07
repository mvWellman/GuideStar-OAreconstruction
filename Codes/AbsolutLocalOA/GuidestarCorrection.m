function [LHS,RHS] = GuidestarCorrection(M1,M2)
% GuidestarCorrection Computes the LHS and RHS correction matrices for
% absolute optic axis measurements using birefringent catheter sheath optic axis orientation
%
%   This function is used for correcting the cumulative round-trip signals
%   from the inside (R1) and outside (R2) of a birefringent sheath, which
%   is assumed to be a Q-oriented linear retarder. The function derives the
%   LHS and RHS multipliers needed to obtain absolute optic axis
%   measurements from these signals. Input:
%       R1 - Cumulative round-trip signal from the inside of the
%       birefringent sheath. R2 - Cumulative round-trip signal from the
%       outside of the birefringent sheath.
%
%   Output:
%       LHS - The left-hand side correction matrix. RHS - The right-hand
%       side correction matrix.

% find Lsys
m1 = decomposeRot(M1);
fun = @(x) findLsys(x,M1);
start = double(mean(m1,2)/2);
xopt = fminsearch(fun,start([1;2]));
[~,m1C] = fun(xopt);

Lsys = makeRot3x3(cat(1,-xopt,0));
inverseRootSquareM1 = makeRot3x3(m1C);

phi = unwrap(atan2(m1C(2,:),m1C(1,:)));

if mean(diff(phi))<0% if negative rotation
    Lsys = makeRot3x3(cat(1,(-xopt+pi*xopt/norm(xopt)),0));
    inverseRootSquareM1 = pagemtimes(Lsys,pagemtimes(M1,Lsys));
end


% correct M2 with the same Lsys
M2C = pagemtimes(Lsys,pagemtimes(M2,Lsys));
m2C = decomposeRot(M2C);
inverseRootSquareM2 = makeRot3x3(-m2C/2);

%take sheath signal
isolatedSheath = pagemtimes(inverseRootSquareM2,pagemtimes(inverseRootSquareM1,inverseRootSquareM2));
SheathSignal = -(unwrapOA1D(isolatedSheath));

% also compute the 2pi wrapped version
SheathSignal2PiWrapped = SheathSignal - 2*pi*SheathSignal./sqrt(sum(SheathSignal.^2,1));

%pick the one with the smaller ret
if mean(sqrt(sum(SheathSignal.^2,1)))>mean(sqrt(sum(SheathSignal2PiWrapped.^2,1)))
    SheathSignal = SheathSignal2PiWrapped;
end

Vtotal = makeRot3x3(cat(1,zeros(2,size(M1,3)),-atan2(SheathSignal(2,:),SheathSignal(1,:))));


RHS = pagemtimes(Lsys,pagemtimes(inverseRootSquareM2,pagetranspose(Vtotal)));
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

