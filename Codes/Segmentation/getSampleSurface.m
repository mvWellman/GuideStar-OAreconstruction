function sampleSurfacePosition = getSampleSurface(inputImage, OuterSheathPosition, intTh, filtw)  
% getSampleSurface - Calculate the sample surface position based on the
% given inputImage matrix, OuterSheathPosition and thresholding.
%
% SYNTAX:
%   sampleSurfacePosition = getSampleSurface(Int, OuterSheathPosition, intTh, filtw)
%
% DESCRIPTION:
%   This function processes an inputImage matrix to determine the position
%   of the sample surface. The function identifies regions in the matrix
%   where the intensity exceeds a specified threshold, and then calculates
%   the position of the surface based on the longest contiguous segment
%   of thresholded pixels. The result is filtered using a median filter
%   to smooth out the positions.
%
% INPUTS:
%   Int - A 2D matrix of intensity values.
%   OuterSheathPosition - Scalar value representing the position of the outer
%                         sheath which acts as a lower bound for surface detection.
%   intTh - Scalar threshold value used to identify significant intensity regions.
%   filtw - Scalar value specifying the width of the median filter window for smoothing.
%
% OUTPUT:
%   sampleSurfacePosition - A vector containing the filtered positions of the sample surface.


% Get the number of rows in the intensity matrix
[dim1, ~] = size(inputImage);

% Create a binary mask where intensity values exceed the threshold
surfaceMask = (inputImage) > intTh;

% Apply the outer sheath position to filter out regions below it
surfaceMask = surfaceMask & (1:dim1)' > OuterSheathPosition;

% Identify the locations where the mask transitions from true to false
edgeDown = circshift(surfaceMask, 1) & ~surfaceMask; % Detect end of thresholded regions

% Compute the length of contiguous segments of thresholded pixels
% diff() calculates the difference between cumulative sums to find segment lengths
% cummax() is used to ensure cumulative maxima for continuity
% max() retrieves the maximum length of contiguous segments
[mv, mp] = max(diff(cummax(cumsum(surfaceMask, 1) .* edgeDown, 1), 1, 1), [], 1);

% Calculate surface positions based on segment lengths
surf = mp - mv;

% Apply median filtering to smooth the surface positions
% Concatenate segments to handle edge cases and apply the median filter
surf = medfilt1(cat(2, surf(end-filtw+1:end), surf, surf(1:filtw)), filtw);

% Output the smoothed surface positions
sampleSurfacePosition = surf(filtw+1:end-filtw);
end

