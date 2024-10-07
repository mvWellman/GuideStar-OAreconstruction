function out = PSProcessGlobalSymmetricCumulative(S1,S2,procStruct)
%out = PSProcessGlobalSymmetricCumulative(S1,S2,procStruct) - computes the 
%cumulative SO3 rotation matrices defined by the Stokes vectors S1 and S2, 
%corresponding to the Stokes vectors measured for an input polarization 
%state modulated between states orthogonal on the Poincaree sphere, using 
%spectral binning. It makes the SO3 matrices D-transpose symmetric and
%aligns them to the central bin, according the systemCopmensation, wich is
%estimated based on the data if not provided by procStruct.

% INPUTS:
%   S1        - 5D matrix of Stokes parameters for channel 1, dimension [Nz, NAlines, NSlices, Nbins, 4].
%   S2        - 5D matrix of Stokes parameters for channel 2, dimension [Nz, NAlines, NSlices, Nbins, 4].
%   procStruct - Structure containing processing parameters:
%                - axial: Vector of axial positions.
%                - scalingFactor: Scaling factor for spectral binning.
%                - fwz: Optional filtering width in z-dimension.
%                - dopTh: Optional threshold for degree of polarization.
%                - roiz: Optional axial range for processing.
%                - systemCompensation: Optional system compensation matrix.
%
% OUTPUTS:
%   out       - Structure containing processed results:
%                - systemCompensation: Compensated system matrix.
%                - dop: Degree of polarization.
%                - S1f: Filtered Stokes parameters for channel 1.
%                - S2f: Filtered Stokes parameters for channel 2.
%                - MM: Measurement matrix of rotation matrices.
%                - MMcorr: (Optional) Corrected measurement matrix.
%

% Define axial range of interest based on scaling factor
locroiz = (floor((procStruct.axial(1)-1)*procStruct.scalingFactor)+1):min(ceil((procStruct.axial(end)-1)*procStruct.scalingFactor)+1,size(S1,1));
S1=S1(locroiz,:,:,:);
S2=S2(locroiz,:,:,:);

dopTh = 0.7;
roiz = 1:size(S1,1);
sysComp=systemCompensation;
% Extract optional parameters from procStruct
fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
    elseif strcmp(fnames{ind},'dopTh')
        dopTh = procStruct.dopTh;
    elseif strcmp(fnames{ind},'roiz')
        roiz = procStruct.roiz;
    elseif strcmp(fnames{ind},'systemCompensation')
        sysComp = procStruct.systemCompensation;
    end
end

% stokesFiltering 
[S1f, S2f, dop, QUVf1,QUVf2] = stokesFiltering(S1,S2, procStruct);

% Normalize the Stokes parameters
S1n = S1f(:,:,:,2:4)./QUVf1;
S2n = S2f(:,:,:,2:4)./QUVf2;


% construct orthonormal tripod for these data points
na = S1n + S2n;
nb = S1n - S2n;
na = na./sqrt(sum(na.^2,4));
nb = nb./sqrt(sum(nb.^2,4));
S1n = (na + nb)/sqrt(2);
S2n = (na - nb)/sqrt(2);

% construct 3x3 rotation matrix at each pixel
S3n = cross(S1n,S2n,4);

% Generate measurement matrix (3x3 rotation matrix) with dimensions [3x3, Nz, NAlines, Nbins]
MM = permute(cat(5,S1n,S2n,S3n),[4,5,1,2,3]);


% Initialize and apply degree of polarization threshold
mm = zeros(size(dop));
mm(roiz, :) = 1;
dopm = dop .* mm;

% Symmetrize matrices and align with respect to the central bin
[MM, sysComp] = compensateSystem(MM, dopm > dopTh,sysComp);
out.systemCompensation = sysComp;

% Store results in the output structure
out.dop = dop;
out.S1f = S1f;
out.S2f = S2f;
out.MM = MM;
% out.MMcorr = MMcorr; % Uncomment if MMcorr is computed

% Interpolate degree of polarization
out.dop = interp1((locroiz - 1)', gather(out.dop), (procStruct.axial - 1)' * procStruct.scalingFactor, 'spline', 0);

% Interpolate the SO3 cumulative round-trip matrix
M = permute(interp1((locroiz - 1)', permute(real(gather(out.MM)), [3, 4, 1, 2]), ...
    (procStruct.axial - 1)' * procStruct.scalingFactor + 1/2, 'spline', 0), [3, 4, 1, 2]);

% Convert the matrix to Euclidean rotation matrices
out.MM = projectToSO3DSymmetric(M);
end




