function strip = getStrip(intensity, firstALine, numberOfALines)
% getStrip - Extract a strip from an intensity image.
%
% SYNTAX:
%   strip = getStrip(intensity, firstALine, numberOfALines)
%
% DESCRIPTION:
%   This function extracts a portion (strip) from a given intensity image based on the specified
%   starting A-line and the number of A-lines to include in the strip. A strip is essentially a
%   collection of contiguous A-lines from the image. This function helps in dividing large images
%   into manageable strips for further processing.

% Calculate the index of the last A-line in the strip.
lastALine = firstALine + numberOfALines - 1;

% Extract the portion of the image corresponding to the strip.
strip.image = intensity(:, firstALine:lastALine);

% Store metadata about the strip.
strip.firstALine = firstALine;
strip.lastALine = lastALine;
strip.intensity = intensity;
end