function [weights,node1,node2] = getNodesAndWeights(paddedTile)
% getNodesAndWeights computes weights and node indices for each edge
% in a padded image tile, based on the normalized gradient image.
%
% INPUTS:
%   paddedTile - 2D matrix representing the image (paddedTile).
%
% OUTPUTS:
%   weights - Vector of calculated weights for edges between nodes.
%   node1   - Vector of indices for the first node in each edge.
%   node2   - Vector of indices for the second node in each edge.
    % Constants for weight calculation
minWeight = 1E-3;   % Minimum weight value
maxWeight = 1000;   % Maximum weight value
% Size of the padded tile
sizePaddedTile = size(paddedTile);

% Initialize matrices for weights and node indices
weight = nan([sizePaddedTile, 8]); % Weights for each of the 8 possible neighbors
node1 = nan([sizePaddedTile, 8]);  % Indices for the first node
node2 = nan([sizePaddedTile, 8]);  % Indices for the second node

% Relative coordinates of the 8 neighboring nodes
neighborX = [1 1 1 0 0 -1 -1 -1];
neighborY = [1 0 -1 1 -1 1 0 -1];


% Compute the normalized gradient of the padded tile.
% The choice between 'negative' or 'positive' depends on the boundary
% characteristics you aim to segment.
normalizedGradientOfPaddedTile = getNormalizedGradImage(paddedTile,'negative');

% Loop through each pixel in the padded tile
for i = 1:sizePaddedTile(1)
    for j = 1:sizePaddedTile(2)
        for z = 1:8
            % Coordinates of the neighboring node
            ii = i + neighborX(z);
            jj = j + neighborY(z);

            % Check if the neighboring node is within the image bounds
            if ii >= 1 && ii <= sizePaddedTile(1) && jj >= 1 && jj <= sizePaddedTile(2)
                % Calculate the weight based on normalized gradient values
                if jj == 1 || jj == sizePaddedTile(2)
                    weight(i, j, z) = minWeight;  % Assign minimum weight for boundary pixels
                elseif normalizedGradientOfPaddedTile(i, j) == 0
                    weight(i, j, z) = maxWeight;  % Assign maximum weight for pixels with zero gradient
                else
                    % Compute weight based on the gradient values of the pixel and its neighbor
                    weight(i, j, z) = 2 - normalizedGradientOfPaddedTile(i, j) - normalizedGradientOfPaddedTile(ii, jj) + minWeight;
                end

                % Calculate linear indices for the nodes
                node1(i, j, z) = sub2ind(sizePaddedTile, i, j);
                node2(i, j, z) = sub2ind(sizePaddedTile, ii, jj);
            end
        end
    end
end
% Reshape matrices to vectors for output
weight = reshape(weight, [sizePaddedTile(1) * sizePaddedTile(2) * 8, 1]);
node1 = reshape(node1, [sizePaddedTile(1) * sizePaddedTile(2) * 8, 1]);
node2 = reshape(node2, [sizePaddedTile(1) * sizePaddedTile(2) * 8, 1]);

% Create a mask to exclude NaN values from the output
mask = ~isnan(weight) & ~isnan(node1) & ~isnan(node2);

% Apply mask to get the final output values
weights = weight(mask);
node1 = node1(mask);
node2 = node2(mask);
end

