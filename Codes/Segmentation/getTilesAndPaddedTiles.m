function [tile,paddedTile]= getTilesAndPaddedTiles(inputImage,tilenumberOfALines,i)
% getTilesAndPaddedTiles - Retrieve a specific tile and its padded version from the intensity image.
%
% SYNTAX:
%   [tile, paddedTile] = getTilesAndPaddedTiles(inputImage, tilenumberOfALines, i)
%
% DESCRIPTION:
%   This function extracts a tile from a given intensity image and pads it with
%   additional columns of zeros. The tile is selected based on its position
%   specified by the index `i`. The function also generates a padded version of the tile
%   to account for edge effects or other processing requirements.

% Determine the total number of A-lines in the input image.
[~, numberofAlines] = size(inputImage);

% Calculate the total number of strips by dividing the number of A-lines by
% the number of A-lines per tile.
numberofstrips = numberofAlines / tilenumberOfALines;

% Retrieve the strips from the intensity image.
strips = getStrips(inputImage, numberofstrips);

% Extract the specific tile from the i-th strip.
tile = getTiles(strips{i});

% Get the size of the extracted tile.
sizeTile = size(tile.image);

% Initialize the padded tile with zeros. The padding is done by adding
% two extra columns (one on each side).
paddedTile = zeros([sizeTile(1), sizeTile(2) + 2]);

% Copy the tile image into the padded tile, centering it in the middle
% columns of the paddedTile matrix.
paddedTile(:, 2:sizeTile(2) + 1) = tile.image;
end



