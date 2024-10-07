function tile = getTileFromStrip(strip, zCenter, dz)
% getTileFromStrip - Extract a tile from a strip based on the specified center and extent.
%
% SYNTAX:
%   tile = getTileFromStrip(strip, zCenter, dz)
%
% DESCRIPTION:
%   This function extracts a tile from a given strip based on the center index and
%   the extent (half-height) of the tile in the Z-dimension. The function adjusts
%   the tile boundaries according to the provided center and extent, and returns
%   the extracted tile along with its metadata.

% Extract the image data for the tile from the strip based on the specified center and extent.
% The tile is defined as a segment centered at "zCenter" with a vertical extent of "dz".
tile.image = strip.image(zCenter - dz : zCenter + dz, :);

% Set the metadata for the tile.
tile.firstZ = zCenter - dz; % The starting index in the Z-dimension for the tile.
tile.lastZ = zCenter + dz;  % The ending index in the Z-dimension for the tile.
tile.zCenter = zCenter;     % The center index in the Z-dimension of the tile.
tile.dz = dz;               % Half the size in the Z-dimension of the tile.
tile.strip = strip;         % The original strip structure from which the tile was extracted.
end

