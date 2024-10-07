function coords = stripToIntensityCoords(strip, stripCoords)
% stripToIntensityCoords - Convert strip coordinates to intensity image coordinates.
%
% SYNTAX:
%   coords = stripToIntensityCoords(strip, stripCoords)
%
% DESCRIPTION:
%   This function converts coordinates from a strip context to the corresponding
%   intensity image coordinates. It adjusts the column index based on the offset
%   provided by the strip data to align with the full intensity image.

% Convert strip coordinates to intensity image coordinates by adjusting
% the column index with the offset provided by the strip structure.
coords = [stripCoords(1), stripCoords(2) + (strip.strip.firstALine-1)];
end



