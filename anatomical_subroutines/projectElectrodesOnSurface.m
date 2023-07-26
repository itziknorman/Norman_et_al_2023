function projectedElectrodes=projectElectrodesOnSurface(electrodes,surfaces,transform)
% project mgrid electrode coordinates on freesurfer/suma surface

space='dicomXYZ';

% projection to the surface

nearestIndex=nan(electrodes.nElec,numel(surfaces));
nearestValues=nan(electrodes.nElec,3,numel(surfaces));
nearestDistance=nan(electrodes.nElec,numel(surfaces));
for iSurface=1:numel(surfaces)
    [nearestIndex(:,iSurface),nearestValues(:,:,iSurface)] = mesh_vertex_nearest(surfaces(iSurface).vertices.(space),electrodes.coord.(space));
    nearestDistance(:,iSurface)=sqrt(sum((nearestValues(:,:,iSurface)-electrodes.coord.(space)).^2,2));
end
[distanceInMM,bestSurfaceInd]=min(nearestDistance,[],2);

bestNearestValues=nan(electrodes.nElec,3);
bestnearestIndex=nan(electrodes.nElec,1);
for iElec=1:electrodes.nElec
    bestNearestValues(iElec,:)=nearestValues(iElec,:,bestSurfaceInd(iElec));
    bestnearestIndex(iElec)=nearestIndex(iElec,bestSurfaceInd(iElec));
end

projectedElectrodes=electrodes;
projectedElectrodes.coord=struct;
projectedElectrodes.coord.(space)=bestNearestValues;
projectedElectrodes.hemisphere={surfaces(bestSurfaceInd).hemisphere}';
projectedElectrodes.nodeInd=bestnearestIndex;
projectedElectrodes.distanceInMMToMesh=distanceInMM;
if isfield(transform,'IJK2XYZ')
    projectedElectrodes.coord.IJK=transform.IJK2XYZ.transformPointsInverse(projectedElectrodes.coord.dicomXYZ);
end

if isfield(transform,'tkregXYZ')
    projectedElectrodes.coord.tkregXYZ=transform.IJK2XYZ_tkreg.transformPointsForward(projectedElectrodes.coord.IJK);
end

if isfield(transform,'XYZ2AFNI')
    projectedElectrodes.coord.afniXYZ=transform.XYZ2AFNI.transformPointsForward(projectedElectrodes.coord.dicomXYZ);
end

end

