function normalizedGradImage = getNormalizedGradImage(Image, varargin)  
% getNormalizedGradImage - Compute and normalize the gradient of an image.
%
% SYNTAX:
%   normalizedGradImage = getNormalizedGradImage(Image, varargin)
%
% DESCRIPTION:
%   This function calculates the gradient of an image and normalizes it. The gradient highlights
%   changes in intensity between adjacent pixels, and normalization scales the gradient values 
%   to a [0, 1] range. The type of normalization can be specified to either emphasize increasing 
%   intensity (positive gradient) or decreasing intensity (negative gradient).

% Get the size of the input image.
sizeImage = size(Image); 

% Initialize the gradient image matrix with NaNs.
gradImg = nan(sizeImage);

% Compute the gradient for each column of the image.
for i = 1:sizeImage(2)
    gradImg(:,i) = gradient(Image(:,i), 2);
end

% Normalize the gradient image based on the specified type.
if nargin > 1 && strcmp(varargin{1}, 'negative')
    % Normalize for negative gradient: emphasizes decreasing intensity.
    normalizedGradImage = (gradImg - min(gradImg(:))) / (max(gradImg(:)) - min(gradImg(:)));
    normalizedGradImage = -1 * normalizedGradImage + 1;
else
    % Normalize for positive gradient: emphasizes increasing intensity.
    normalizedGradImage = (gradImg - min(gradImg(:))) / (max(gradImg(:)) - min(gradImg(:)));
end
end
