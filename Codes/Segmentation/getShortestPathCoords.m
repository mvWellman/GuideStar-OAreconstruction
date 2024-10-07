function [pathCoordX,pathCoordY] = getShortestPathCoords(G,Tile)
% getShortestPathCoords - Compute the shortest path coordinates in a tile
% based on a graph.
%
% SYNTAX:
%   [pathCoordX, pathCoordY] = getShortestPathCoords(G, Tile)
%
% DESCRIPTION:
%   This function finds the shortest path through a graph representation of
%   a padded tile image. The function computes the shortest path from the
%   top-left to the bottom-right of the padded image and converts the
%   resulting path indices into coordinates. It also ensures the path
%   contains only vertical changes by filtering out horizontal transitions.

% Get the size of the tile image.
sizeTile = size(Tile.image);

% Pad the tile image with two columns of zeros on the left and right.
paddedTile = zeros([sizeTile(1), sizeTile(2) + 2]);
paddedTile(:, 2:sizeTile(2) + 1) = Tile.image;

% Compute the shortest path from the top-left to the bottom-right of the
% padded image.
pathInd = shortestpath(G, 1, numel(paddedTile(:)));

% Convert the linear indices of the path to subscript coordinates.
[pathCoordX, pathCoordY] = ind2sub(size(paddedTile), pathInd);

% Filter out the coordinates to include only vertical changes in the path.
pathCoordX = pathCoordX(gradient(pathCoordY) ~= 0);
pathCoordY = pathCoordY(gradient(pathCoordY) ~= 0);

% Adjust the coordinates to match the original tile size (excluding
% padding).
pathCoordX = pathCoordX(1:sizeTile(2));
pathCoordY = pathCoordY(1:sizeTile(2));
end
