function r = unwrapOA1D(Rin)
% unwrap rotation vector of SO3 linear retarder Rin along the third dimension.
% The first two dimensions of Rin have to be 3x3

dim = size(Rin);
dim = dim(2:end);

% eliminate additional dimensions, if present
Rin = Rin(:,:,:,:);
% get rotation vector
r = decomposeRot(Rin);
Rsqinv = makeRot3x3(-r/2);

% generate zeros to concatenate
dd = dim;
dd(2) = 1;
zz = zeros(dd);
mx = cat(2,zz,decomposeRot(pagemtimes(Rsqinv(:,:,1:end-1,:),pagemtimes(Rin(:,:,2:end,:),Rsqinv(:,:,1:end-1,:)))));
dd = cat(2,zz,r(:,2:end,:)-r(:,1:end-1,:));

cp =(sum(mx.*dd,1)./sqrt(sum(mx.^2,1))./sqrt(sum(dd.^2,1)));
mm = mod(cumsum(cp<0,2),2)>0;

% and its 2*pi wrapped version
rw = r - 2*pi*r./sqrt(sum(r.^2,1));

r(:,mm) = rw(:,mm);

r = reshape(r,dim);
