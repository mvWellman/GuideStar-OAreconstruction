function depthResolvedRotationMatrix = getDepthResolvedRet(CumulativeRotationMatrix,structure)
% getDepthResolvedRet Computes depth-resolved retardance or rotation from a
% cumulative retardance matrix.
%
% depthResolvedRotationMatrix =
% getDepthResolvedRet(CumulativeRotationMatrix, structure) takes a
% cumulative rotation matrix and a structure containing various parameters
% to compute the depth-resolved rotation matrix.
%
% Inputs:
%   CumulativeRotationMatrix - The cumulative retardance matrix.. structure
%   - A structure containing the following fields:
%     - dop: Degree of polarization (DOP) values. - dopTh: Threshold for
%     DOP. - surf: Surface positions or depths where the sample surface
%     starts.
%
% Output:
%   depthResolvedRotationMatrix - Depth-resolved rotation matrices.

fnames = fieldnames(structure);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'dop')
        dop = structure.dop;
    elseif strcmp(fnames{ind},'dopTh')
        dopTh = structure.dopTh;
    elseif strcmp(fnames{ind},'surf')
        surf = structure.surf;
    end
end

% Create a mask to filter out data below the surface position and where DOP
% is below threshold
mask = dop>dopTh;
mask = mask & (1:size(mask,1))'>surf;

% Initialize the depth-resolved rotation matrix with identity matrices.
% N is a 4D matrix with dimensions [3, 3, 1, depth], where depth is the 4th
% dimension of CumulativeRotationMatrix
N = repmat(eye(3),[1,1,1,size(CumulativeRotationMatrix,4)]);

% Iterate over each depth position to compute the depth-resolved retardance
    for indz = 1:size(CumulativeRotationMatrix,3)
        
        % Compute local rotation matrix nloc by applying the current correction matrix N
        Nloc = pagemtimes(N.*[1,1,-1;1,1,-1;-1,-1,1],pagemtimes(CumulativeRotationMatrix(:,:,indz,:),pagetranspose(N)));
        % Decompose the local rotation matrix to extract the rotation components
        Nloc = decomposeRot(Nloc);
        locRotMatrix(:,indz,:) = Nloc;  % Store the decomposed rotation components
        locRotMatrix(:,indz,~mask(indz,:)) = 0; %Zero out where mask is false

        % Compute the rotation matrix for the next depth
        NlocRoot = makeRot3x3(Nloc/2);
        Nnext = pagemtimes(NlocRoot,N);
        % Update the identity matrix with the next depth's rotation matrix where not masked
        N(:,:,:,mask(indz,:)>0) = Nnext(:,:,:,mask(indz,:)>0);% only multipy with Nnext when not masked
        
    end
    % Return the depth-resolved rotation matrix
    depthResolvedRotationMatrix=locRotMatrix;
end