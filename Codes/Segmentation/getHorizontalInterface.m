function horizontalPosition = getHorizontalInterface(Tile)
% getHorizontalInterface - Calculate the horizontal positions of the
% interface within a given tile.
%
% SYNTAX:
%   horizontalPosition = getHorizontalInterface(Tile)
%
% DESCRIPTION:
%   This function processes a tile of an image to determine the horizontal
%   positions of the interface. It first pads the tile image with two
%   columns of zeros on the sides, computes the gradient image, finds nodes and
%   weights for a graph reconstruction, and then calculates the shortest
%   path to determine interface positions.

% INPUTS:
%   Tile - A structure containing all information needed to recover the
%   position of a tile in the main inputImage, e.g.:
%     - image: A 2D matrix representing the image data of the tile.
%
% OUTPUTS:
%   horizontalPosition - A matrix where each row represents the (x, y)
%   coordinates of the interface positions in the tile.

% Get the size of the tile image
sizeTile = size(Tile.image);

% Pad the tile image with zeros on the left and right sides
paddedTile = zeros([sizeTile(1), sizeTile(2) + 2]);
paddedTile(:, 2:sizeTile(2) + 1) = Tile.image;

% Compute nodes and weights for graph reconstruction
[weights, nodes1, nodes2] = getNodesAndWeights(paddedTile);

% Create a graph from the nodes and weights
G = graph(nodes1(:), nodes2(:), weights(:));

% Compute the shortest path coordinates using the graph
[pathCoordX, pathCoordY] = getShortestPathCoords(G, Tile);

% Initialize the output matrix for horizontal positions
horizontalPosition = [];

% Convert path coordinates from the padded tile to original inputImage coordinates
for n = 1:numel(pathCoordX)
    % Convert (x, y) coordinates from padded tile to original inputImage
    coords = tileToIntensityCoords(Tile, [pathCoordX(n), pathCoordY(n)]);
    % Append the coordinates to the result matrix
    horizontalPosition = cat(1, horizontalPosition, coords);
end

end

