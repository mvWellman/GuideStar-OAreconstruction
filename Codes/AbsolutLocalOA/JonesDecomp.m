function [r,d] = JonesDecomp(J,polar,pureDiatt)
% JonesDecomp computes the retardation and diattenuation from a Jones
% matrix, provided as input argument ( the matrix logarithm is part of the
% function, using the 'concurrent' decomposition.
% If called with the additioanl argument 'polar', it performs the polar
% decomposition, decomposing the Jones matrix into a sequence of a
% diattenuation and a retardation matrix, J = Jr*Jd.
% Different approach without logm
% J is either 2x2xN, where the first 4 elements are [J11;J21;J12;J22]

dim = size(J);

if numel(dim)<=2
    dim(3) = 1;
end

if nargin<2
    polar = false;
end
if nargin<3
    pureDiatt = false;
end

if ~polar
    detJ = sqrt(J(1,1,:).*J(2,2,:)-J(2,1,:).*J(1,2,:));
    J = J(:,:,:)./detJ;
    q = cat(1,(J(1,1,:)-J(2,2,:)),(J(2,1,:)+J(1,2,:)),(-1i*J(2,1,:)+1i*J(1,2,:)))/2;
    tr = (J(1,1,:) + J(2,2,:))/2;
    c = acosh(tr);
    csin = c./sinh(c);
    csin(c==0) = 1;
    f = 2*q.*csin;
    r = reshape(-imag(f),[3,dim(3:end)]);
    d = reshape(real(f),[3,dim(3:end)]);
else% polar decomposition
    detJ = sqrt(J(1,1,:).*J(2,2,:)-J(2,1,:).*J(1,2,:));
    J = J(:,:,:)./detJ;
    [~,d] = JonesDecomp(pagemtimes(pagectranspose(J),J));

    d = d/2;
    r = JonesDecomp(pagemtimes(J,makeJones(zeros(size(d)),-d)));
    r = reshape(r,[3,dim(3:end)]);
    d = reshape(d,[3,dim(3:end)]);
end
