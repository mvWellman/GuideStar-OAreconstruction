function coords = tileToIntensityCoords(tile, tileCoords)
% tileToIntensityCoords - Convert tile coordinates to intensity coordinates.
%
% SYNTAX:
%   coords = tileToIntensityCoords(tile, tileCoords)
%
% DESCRIPTION:
%   This function converts coordinates from a given tile to the corresponding
%   intensity coordinates in the full image. The conversion is performed in two
%   steps: first by converting tile coordinates to strip coordinates, and then
%   by converting those strip coordinates to intensity coordinates.


% Convert tile coordinates to strip coordinates.
stripCoords = tileToStripCoords(tile, tileCoords);

% Convert strip coordinates to intensity coordinates.
coords = stripToIntensityCoords(tile, stripCoords);

end