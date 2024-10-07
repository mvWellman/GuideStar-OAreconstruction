function strips = getStrips(intensity, numberOfStrips)
% getStrips - Divide an intensity image into a specified number of strips.
%
% SYNTAX:
%   strips = getStrips(intensity, numberOfStrips)
%
% DESCRIPTION:
%   This function divides an intensity image into a specified number of strips. Each strip is a
%   collection of A-lines extracted from the input image. The function handles the case where
%   the number of A-lines is not perfectly divisible by the number of strips by distributing
%   the leftover A-lines across the strips.

% Get the number of A-lines in the intensity image.
[~, numberOfALines] = size(intensity);

% Calculate the number of A-lines per strip and handle any leftover A-lines.
stripSize = floor(numberOfALines / numberOfStrips);
leftoverALines = mod(numberOfALines, numberOfStrips);

% Initialize the sizes for each strip.
if leftoverALines == 0
    stripSizes = repmat(stripSize, 1, numberOfStrips);
else
    stripSizes = repmat(stripSize, 1, numberOfStrips - 1);
    stripSizes(end + 1) = stripSize + leftoverALines;
end

% Initialize the cell array to hold the strips.
strips = cell(1, numberOfStrips);

% Index for the starting A-line of the current strip.
start_of_strip = 1;

% Loop to create each strip from the intensity image.
for k = 1:numberOfStrips
    % Extract the strip from the intensity image.
    strip = getStrip(intensity, start_of_strip, stripSizes(k));
    
    % Store the extracted strip in the cell array.
    strips{k} = strip;
    
    % Update the starting A-line index for the next strip.
    start_of_strip = start_of_strip + stripSizes(k);
end
end

