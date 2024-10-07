function tile = getTiles(strip)
% getTiles - Extract tiles from strip and prepare its metadata.
%
% SYNTAX:
%   tile = getTiles(strip)
%
% DESCRIPTION:
%   This function extracts a tile from a strip and prepares its metadata.
%   The tile is essentially a segment of the strip with additional metadata
%   about its size and position within the strip. This information is used
%   for processing and further analysis of the tile within the larger image context.
% INPUTS:
%   strip - A structure containing information about the strip, including:
%     - `image`: A 2D matrix representing the image data for the strip.
%
% OUTPUTS:
%   tile - A structure containing the extracted tile's image and associated metadata, including:
%     - `image`: The image data of the tile.
%     - `zCenter`: The center index of the Z-dimension of the tile.
%     - `dz`: Half of the size in the Z-dimension of the tile, used to define the extent of the tile.
%     - `firstZ`: The starting index of the Z-dimension for the tile.
%     - `lastZ`: The ending index of the Z-dimension for the tile.
%     - `strip`: The original strip structure from which the tile was extracted.

% Retrieve the image data from the strip structure.
im = strip.image;

% Determine the size of the image.
sizeIm = size(im);

% Initialize the tile structure.
tile.image = im;

% Calculate the center index in the Z-dimension (vertical axis) of the tile.
tile.zCenter = (sizeIm(1) + 1) / 2;

% Calculate half of the size in the Z-dimension, which is used to define the extent of the tile.
tile.dz = (sizeIm(1) - 1) / 2;

% Define the starting and ending indices in the Z-dimension for the tile.
tile.firstZ = 1;
tile.lastZ = sizeIm(1);

% Include the original strip information in the tile structure.
tile.strip = strip;
end

