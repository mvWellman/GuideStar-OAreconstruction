function [innerSheathPosition, outerSheathPosition, sampleSurfacePosition] = getSheathAndSampleInterfacePosition(inputImage,innerSheath,outerSheath,varargin)

% Define default values for optional parameters
p = inputParser;

addParameter(p, 'dz', 20);                 % Default for dz
addParameter(p, 'dzMaskOuterSheath', 5);    % Default for ZmaskOuterSheath
addParameter(p, 'tileNumberOfALines', 32); % Default for tileNumberOfALines

% Parse input arguments
parse(p, varargin{:});

% Assign parsed values to variables
dz = p.Results.dz;
dzMaskOuterSheath = p.Results.dzMaskOuterSheath;
tileNumberOfALines = p.Results.tileNumberOfALines;


% getSheathAndSampleInterfacePosition - Determine the positions of the
% sheath and sample interfaces in catheter-based OCT images.
%
% SYNTAX:
%   [innerSheathPosition, outerSheathPosition, sampleSurfacePosition] =
%   getSheathAndSampleInterfacePosition(inputImage, innerSheath,
%   outerSheath, tileNumberOfALines)
%
% DESCRIPTION:
%   This function processes an input image to identify and return the
%   positions of the inner sheath, outer sheath, and sample surface
%   interfaces using Graph Cut Segmentation (GTS). In the OCT paradigm, the
%   image dimensions are organized such that the second dimension
%   (horizontal axis) represents the number of A-lines, while the first
%   dimension (vertical axis) denotes the number of pixels in each A-line.
%   The positions are calculated based on the provided inputs.

% INPUTS
%   inputImage - A 2D numeric array representing the image from which the
%   positions will be extracted. innerSheath - Inner sheath position in the
%   first A-line of entire image. outerSheath - Outer sheath position in
%   the first A-line of entire image. tileNumberOfALines - The
%   tileNumberOfALines parameter is used to divide the entire OCT image
%   into smaller sub-images, referred to as Tiles. This parameter specifies
%   the number of A-lines contained within each Tile. dz - Axial range to
%   limit the search area around the inner sheath. ZmaskOuterSheath - A mask
%   range for outer sheath detection

% OUTPUTS:
%   innerSheathPosition - A vector containing the positions of the inner
%   sheath interface in the image. outerSheathPosition - A vector
%   containing the positions of the outer sheath interface in the image.
%   sampleSurfacePosition - A vector containing the positions of the sample
%   surface in the image.


% Get the number of A-lines in the image
[~, numberofAlines] = size(inputImage);
% Calculate the number of tiles based on the number of A-lines per tile
numberOfTiles = numberofAlines / tileNumberOfALines;

% Initialization for storing positions
positions = [];
innerSheathPosition = [];

% Define the initial zCenter for the inner sheath and a range around it
zCenter = innerSheath;


% Process each tile to find the inner sheath positions
for i = 1:numberOfTiles
    % Extract the current tile from the entire image
    Tile = getTileFromEntireImage(inputImage, tileNumberOfALines, dz, zCenter, i);
    % Find horizontal interface positions in the current tile
    interfacePosition = getHorizontalInterface(Tile);
    % Append the results to the positions list
    positions = cat(1, positions, interfacePosition(:, 1));
    % Update zCenter for the next tile
    zCenter = positions(end, 1) + 1;
end

% Calculate inner sheath positions by adjusting the positions. The offset
% of -2 shifts the positions from the upper edge of the inner sheath
% interface to its center.
innerSheathPosition = (cat(1, innerSheathPosition, positions) - 2)';

% Initialize variables for outer sheath processing
positions = [];
outerSheathPosition = [];

% Apply median filtering to the input image
filteredIntensity = medfilt2(inputImage, [1 51]);

% Calculate the expected outer sheath positions based on the innerSheath
sheathThickness = outerSheath - innerSheath;
expectedOuterSheathPosition = innerSheathPosition + sheathThickness;

% Create masked image regions around the expected outerSheath positions
for i = 1:numberofAlines
    maskedImage(:, i) = filteredIntensity((expectedOuterSheathPosition(i) - dzMaskOuterSheath):(expectedOuterSheathPosition(i) + dzMaskOuterSheath), i);
end

% Process each tile to find the outer sheath positions
for i = 1:numberOfTiles
    [Tiles, ~] = getTilesAndPaddedTiles(maskedImage, tileNumberOfALines, i);
    interfacePosition = getHorizontalInterface(Tiles);
    positions = cat(1, positions, interfacePosition(:, 1));
end

% Adjust outer sheath positions based on the differences from expected
% positions
positions = positions - 2; % -2 due to the same reason as innerSheath adjustment
diff = (cat(1, outerSheathPosition, positions))';
outerSheathPosition = expectedOuterSheathPosition - dzMaskOuterSheath + diff;

% Define parameters for sample surface detection
intTh = 70; % Intensity threshold for sample surface detection
filtw = 20; % Width of median filter for smoothing

% Compute the sample surface position based on the outerSheath positions
sampleSurfacePosition = getSampleSurface(inputImage, outerSheathPosition, intTh, filtw);
end
