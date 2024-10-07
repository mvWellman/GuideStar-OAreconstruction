function R = rotationMatrixAtPosition(MM,positions)
%   R = rotationMatrixAtPosition(MM, positions) returns a rotation matrix
%   extracted from the 4D matrix MM based on the given positions.
%
%   Input:
%       MM - A 4D matrix of size [3, 3, N, M], where N and M are the dimensions
%            of the matrix in the 3rd and 4th dimensions, respectively. Each
%            3x3 matrix in MM represents a rotation matrix.
%       positions - A vector or array specifying the positions from which to
%                   extract the rotation matrices. The values should be in
%                   the range [1, N*M], where each position corresponds to a
%                   unique 3x3 rotation matrix in MM.
%
%   Output:
%       R - A 3x3xP matrix containing the extracted rotation matrices. P is
%           the number of specified positions.


% Calculate the linear indices based on the positions and the size of MM
lInd = (round(positions')-1)*9+1 + (0:size(MM,4)-1)'*9*size(MM,3);
 % Reshape and adjust indices for extracting 3x3 matrices
lInd = shiftdim(lInd(:),-2) + reshape((0:8)',[3,3]);
% Extract rotation matrices from MM using the computed indices
R = (MM(lInd));
end
