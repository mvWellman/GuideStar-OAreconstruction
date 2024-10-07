function [S1, S2, dop, QUVf1, QUVf2] = stokesFiltering(S1,S2, procStruct)

% [S1, S2, dop, QUVf1, QUVf2] = stokesProcessing(S1,S2, procStruct)
% 
% Inputs: S1 and S2 (stokes vectors), procStruct (includes filtering info)
%   S1/S2: matrix (3 or 4 dimensions)
%       4 dimensions -> nZ x nAlines x nSpectralBins x 3 or 4
%       3 dimensions -> No nSpectralBins -> nZ x nAlines x 3 or 4
%       Final dimension size: if 3 -> [U V W] (no intensity), 
%           if 4 -> [I U V W]
%   procStruct: structure
%       Field 'fwx': scalar, lateral filtering
%       Field 'fwz': scalar, vertical filtering
% 
% Outputs: Filtered stokes vectors and normalized degree of polarization
%   S1/S2: matrix 
%       Dimensions -> nZ x nAlines x nSpectralBins x 4 ([I U V W])
%   QUVf1/QUVf2: matrix
%       Dimensions -> nZ x nAlines x nSpectralBins
%   dop: matrix
%       Dimensions -> nZ x nAlines
% 
% The goal of this function is to: 1) calculate the intensity of the Stokes
% vectors if necessary, 2) apply a Gaussian kernel to filter the Stokes
% vectors, and 3) compute the degree of polarization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SETUP
% % INPUT: procStruct
% % OUTPUT: procStruct

if ~isfield(procStruct,'fwx')
    fwx = [];
else
    fwx = procStruct.fwx;
end
if ~isfield(procStruct,'fwz')
    fwz = [];
else
    fwz = procStruct.fwz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANAGE INPUT AND COMPUTE INTENSITY
% % INPUT: S1, S2
% % OUTPUT: S1, S2

dim = size(S1);
% manage case of no spectral binning
if numel(dim)<4 % no spectral binning used
    S1 = permute(S1,[1,2,4,3]);% introduce third dimension
    S2 = permute(S2,[1,2,4,3]);% introduce third dimension
    dim = size(S1);
end

% Compute intensity
if dim(4) == 3 % 3-component Stokes vector was provided
    S1 = cat(4,sqrt(sum(S1.^2,4)),S1);
    S2 = cat(4,sqrt(sum(S2.^2,4)),S2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN FILTERING
% % INPUT: S1, S2, fwx, fwz
% % OUTPUT: S1, S2

if ~isempty(fwx)
    sigma = fwx / sqrt(8) / sqrt(log(2));
    hX = fspecial('gaussian', [1 , 2*ceil(2*sigma)+1], sigma);
    
    if ~isempty(fwz)
        sigmaZ = fwz / sqrt(8) / sqrt(log(2));
        hZ = fspecial('gaussian', [2*ceil(2*sigmaZ)+1,1], sigmaZ);
        h = hZ * hX;
    else
        h = hX;
    end
    
    h = h / sum(h(:)); 
    S1 = imfilter(S1,h,'circular');
    S2 = imfilter(S2,h,'circular');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZATION OF DOP
% % INPUT: S1, S2
% % OUTPUT: QUVf1, QUVf2, dop1, dop2, dop

QUVf1 = sqrt(sum(S1(:,:,:,2:4).^2,4));
QUVf2 = sqrt(sum(S2(:,:,:,2:4).^2,4));

dop1 = QUVf1./S1(:,:,:,1);
dop2 = QUVf2./S2(:,:,:,1);
dop = 1/2*mean(dop1,3) + 1/2*mean(dop2,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end