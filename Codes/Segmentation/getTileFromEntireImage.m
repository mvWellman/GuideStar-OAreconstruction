function tile= getTileFromEntireImage(inputImage,tilenumberOfALines,dz,zCenter,i)
% getTileFromEntireImage - Extract a specific tile from an entire image based on the provided parameters.
%
% SYNTAX:
%   tile = getTileFromEntireImage(inputImage, tilenumberOfALines, dz, zCenter, i)
%
% DESCRIPTION:
%   This function extracts a tile from a larger image. The image is divided into strips based
%   on the number of A-lines per tile, and each strip is then used to extract the specified tile.
%   The function calculates which strip contains the tile and then extracts the tile based on the
%   center and extent parameters provided.

% Determine the total number of A-lines in the input image.
[~, numberofAlines] = size(inputImage);

% Calculate the number of strips required based on the number of A-lines per tile.
numberOfTiles = numberofAlines / tilenumberOfALines;

% Split the entire image into strips based on the number of A-lines per strip.
strips = getStrips(inputImage, numberOfTiles);

% Extract the tile from the appropriate strip using the provided center and extent parameters.
tile = getTileFromStrip(strips{i}, zCenter, dz);
end


 
