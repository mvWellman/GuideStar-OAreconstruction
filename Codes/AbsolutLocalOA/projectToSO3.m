function Rout = projectToSO3(Min,boolDSym)
% Rout = projectToSO3(Min) takes Min, an array of 3x3 matrices (the first 
% 'page') and 'projects' them into SO3. If boolDSym is set to true, the 
% input is assumed to be D-transpose symmetric, were D = diag([1,1,-1]), 
% and D*Min.'*D = Min. boolDsym defaults to false. If true,the output 
% SO3 matrix is known to be a linear retarder, and the projection avoids 
% the trivial solution of Rout = diag([-1,-1,1]).
%
% Conceptually, the projection of Min into SO3 is similar to the polar
% decomposition: Rout = Min/sqrtm(Min.'*Min). However, taking matrix roots
% is ambiguous, and we operate directly on the singular values of Min to 
% ensure SO3 output.
%
% We have Min = u*d*v', where [u,d,v] = svd(Min). The lsq solution to
% fitting Min with an SO3 matrix leads to max trace(Min'*Rout) = max
% trace(v*d*u'*Rout), from where the solution is Rout = u*v'. However,
% det(Min) can be negative, but we need det(Rout)=1. If Min has negative
% determinant, we need to flip the smallest (last) singular value.
%
% If D*Min.'*D = Min, D*Min is a transpose-symmetric matrix and has the
% eigenvalue decomposition D*Min = v*lam*v.', where lam is real but possibly
% negative valued, and abs(lam) = d = lam*S, where S is a diagonal matrix 
% defining negative or positive signs on the diagonal. We need to avoid
% the combination of S that leads to the trivial solution of Rout = D, 
% which is not a linear retarder. If D results in the smallest lsq error, 
% we have to flip the sign of two smallest (last) eigenvalues.

% default of boolDSym
if nargin<2
    boolDSym = false;
end

% register the dimension of Min to enforce the same dimension on the output
dim = size(Min);

% take svd of input matrix
[u,d,v] = pagesvd(Min,'vector');

% compute the matrix determinant
mdet = sign(Min(1,1,:).*Min(2,2,:).*Min(3,3,:) + Min(2,1,:).*Min(3,2,:).*Min(1,3,:) + Min(3,1,:).*Min(1,2,:).*Min(2,3,:) - Min(3,1,:).*Min(2,2,:).*Min(1,3,:) - Min(1,1,:).*Min(3,2,:).*Min(2,3,:) - Min(2,1,:).*Min(1,2,:).*Min(3,3,:));
dd = [ones(size(mdet));ones(size(mdet));mdet];% Rout =  u*dd*v.', where dd = diag(1,1,det(M)) to ensure det(Rout) = 1

if boolDSym % manage symmetric case
    S = pagemtimes(pagetranspose(u),[1;1;-1].*v);% u*d*v' = D*v*lam*v.' = D*v*S*d*v.' => u = D*v*S; hence we can recover S directly from u and v
    mm = (S(1,1,:)+S(2,2,:)+S(3,3,:).*mdet)<-2;% whenever sign(S).*[1,1,mdet] = [-1,-1,-1], we neet to flip the last two elements of dd. Because of numerical accuracy limitations ==-3 does not work. However, the sum only takes the value -3 or 1, and -2 is a robust cutoff.
    dd(2:3,1,mm) = -dd(2:3,1,mm);
end

Rout = reshape(pagemtimes(u(:,:,:),dd.*pagetranspose(v(:,:,:))),dim);% Rout =  u*dd*v.', where dd = diag(1,1,det(M))
