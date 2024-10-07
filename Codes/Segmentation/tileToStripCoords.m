function coords = tileToStripCoords(tile, tileCoords )
% tileToStripCoords - Convert tile coordinates to strip coordinates.
%
% SYNTAX:
%   coords = tileToStripCoords(tile, tileCoords)
%
% DESCRIPTION:
%   This function converts coordinates from a given tile to the corresponding
%   coordinates in the full strip. The conversion is based on the
%   position of the tile within the larger image or strip, as specified by the
%   `tile` structure.

% Convert local tile coordinates to global coordinates by adding the tile's
% vertical offset (firstZ - 1) to the row index, while the column index remains unchanged.
coords = [(tile.firstZ - 1)+ tileCoords(1), tileCoords(2) ];
end

