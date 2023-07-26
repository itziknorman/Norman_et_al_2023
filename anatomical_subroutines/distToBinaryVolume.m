function [numOfAdjVoxels,distanceToVol]=distToBinaryVolume(elecCoord,binaryVolume,transform,radius)

% elecCoord = native space coordinates of an electrode
% binary volume is a 0/1 matrix in ijk space
% transform is the tranfromation matrix from IJX to freesurfer RAS (tkregister)
nonZeroInd=find(binaryVolume);
% linear indecis to ijk coordinates
nonZeroIJK=nan(numel(nonZeroInd),3);
[nonZeroIJK(:,1),nonZeroIJK(:,2),nonZeroIJK(:,3)]=ind2sub(size(binaryVolume),nonZeroInd);
% convert ijk to xyz
nonZeroDicomXYZ = transform.IJK2XYZ_tkreg.transformPointsForward(nonZeroIJK);

if ~exist('radius','var'), radius = 3; end

distanceToVol = nan(1,size(elecCoord,1));
numOfAdjVoxels = nan(1,size(elecCoord,1));

for iElec=1:size(elecCoord,1)
    if size(nonZeroDicomXYZ,1)>0
    % find how many voxels of a particular region are
    % within [radius] from the electrode:
    curCoord = elecCoord(iElec,:); 
    tmp=sqrt(sum(bsxfun(@minus,nonZeroDicomXYZ,curCoord).^2,2));
    nVoxels=sum(tmp<radius);
    numOfAdjVoxels(iElec)=nVoxels;
    % compute min distance to current ROI:
    distanceToVol(iElec)=min(sqrt(sum(bsxfun(@minus,nonZeroDicomXYZ,curCoord).^2,2)));
    else
        numOfAdjVoxels(iElec) = nan;
        distanceToVol(iElec) = nan;
    end
end
    