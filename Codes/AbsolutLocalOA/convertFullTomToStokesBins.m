function [S1,S2,scalingFactor] = convertFullTomToStokesBins(tom,Nmapping,binFract)
%[S1,S2] = convertFullTomToStokesBins(tom) takes a fully coherent
%reconstructed tomogram, containing polarization diverse channels and
%alternated input polarization and converts it into spectral bins of width
%1/binFract. Nmapping is the number of points used in the mapping,
%traditionally half of the number of sampled points.

% SYNTAX:
%   [S1, S2, ssFact] = convertFullTomToStokesBins(tom, Nmapping, binFract)
%
% DESCRIPTION:
%   This function processes a fully coherent reconstructed tomogram, which includes 
%   polarization diverse channels and alternated input polarizations, to convert it 
%   into spectral bins of width 1/binFract. The output is binned according to 
%   the specified fraction of the full spectral width. 
%   The function calculates Stokes parameters for each bin and returns the 
%   results along with the scaling factor for the binning process.
% INPUTS:
%   tom      - A structure containing the fully coherent reconstructed tomogram with 
%              fields `ch1` and `ch2` representing polarization channels.
%   Nmapping - The number of points used in the mapping, traditionally half of the 
%              number of sampled points.
%   binFract  - Fraction representing the width of each spectral bin relative to 
%               the total spectral range. Determines the number of bins.
%
% OUTPUTS:
%   S1       - A 5D matrix of Stokes parameters for channel 1, binned into spectral 
%              bins. The dimensions are (2*W, NAlines, NSlices, Nbins, 4).
%   S2       - A 5D matrix of Stokes parameters for channel 2, binned into spectral 
%              bins. The dimensions are (2*W, NAlines, NSlices, Nbins, 4).
%   scalingFactor   - Scaling factor used in the binning process.


% Get dimensions of the input tomogram
Nz = size(tom.ch1, 1);          % Number of axial points
NAlines = size(tom.ch1, 2) / 2; % Number of A-lines (half since data is alternated)
NSlices = size(tom.ch1, 3);     % Number of slices

% Construct the inverse of the Hanning filter originally applied
whinv = circshift(cat(1, 1 ./ hanning(Nmapping), zeros(Nz - Nmapping, 1)), -Nmapping / 2);

% Calculate start indices for each bin
Nbins = 2 * binFract - 1; % Number of bins
W = round(Nmapping / binFract); % Width of each bin
startIndices = round((0:Nbins-1) * Nmapping / (Nbins + 1));

% Construct the Hanning window for binning
wh = hanning(W);

% Apply inverse filter and shift
sp1 = circshift(ifft(tom.ch1) .* whinv, Nmapping / 2);
sp2 = circshift(ifft(tom.ch2) .* whinv, Nmapping / 2);

% Initialize output matrices
S1 = zeros(2 * W, NAlines, NSlices, Nbins, 4, 'like', sp1);
S2 = zeros(2 * W, NAlines, NSlices, Nbins, 4, 'like', sp2);

% Process each bin
for indw = 1:Nbins
% Apply window and perform FFT
chh1 = fft(sp1(startIndices(indw) + (1:W), 1:2:end, :) .* wh, 2 * W);
chv1 = fft(sp2(startIndices(indw) + (1:W), 1:2:end, :) .* wh, 2 * W);
chh2 = fft(sp1(startIndices(indw) + (1:W), 2:2:end, :) .* wh, 2 * W);
chv2 = fft(sp2(startIndices(indw) + (1:W), 2:2:end, :) .* wh, 2 * W);

% Compute Stokes parameters
S1(:, :, :, indw, :) = cat(5, abs(chh1).^2 + abs(chv1).^2, ... % I
                              abs(chh1).^2 - abs(chv1).^2, ... % Q
                              2 * real(chh1 .* conj(chv1)), ... % U
                              -2 * imag(chh1 .* conj(chv1))); % V
S2(:, :, :, indw, :) = cat(5, abs(chh2).^2 + abs(chv2).^2, ... % I
                              abs(chh2).^2 - abs(chv2).^2, ... % Q
                              2 * real(chh2 .* conj(chv2)), ... % U
                              -2 * imag(chh2 .* conj(chv2))); % V
end

% Squeeze dimensions if only one slice is present
if NSlices == 1
S1 = squeeze(S1);
S2 = squeeze(S2);
end

% Compute the scaling factor for binning
%W: Width of each bin.
%nZ: number of axial points.
scalingFactor = 2 * W / Nz;

end
