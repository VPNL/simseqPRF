function out = simseq_checkVoxelCorrespondance(coords, selectedpRFs)
% Function check nr voxels of BOLD data and pRFs and make sure they correspond

% Check which hemispheres were dealing with
hemi = {};
if isfield(coords.pRFs,'lh')
    hemi{1} = 'lh';
end
if isfield(coords.pRFs,'rh')
    if isempty(hemi), hemi{1} = 'rh';
    else  hemi{2} = 'rh'; end
end

% Extra info from data structs
for h = 1:length(hemi)
    
    currHemi = hemi{h};
    
    % nr of all pRFs, data voxels, pRFs surviving var expl thresholds
    pRFs.totalNrVoxels.(currHemi)  = size(coords.pRFs.(currHemi),2);
    pRFs.inclByVarExpl.(currHemi)  = find(selectedpRFs.predVEmask.(currHemi));
    pRFs.exclByVarExpl.(currHemi)  = find(~selectedpRFs.predVEmask.(currHemi));
    dataRuns.totalNrVoxels.(currHemi) = size(coords.dataRuns.(currHemi),2);
    
    if isequal(find(selectedpRFs.predVEmask.(currHemi)),pRFs.inclByVarExpl.(currHemi))
        coords.pRFs.(currHemi)   = coords.pRFs.(currHemi)(:,pRFs.inclByVarExpl.(currHemi));
    end    
    
    % Check for any differences in voxels (sometimes the pRF is defined,
    % but the voxel has no data (falls outside the FoV) or vice versa.
    [ai, bi, ci] = intersectCols(coords.dataRuns.(currHemi), coords.pRFs.(currHemi));
    dataRuns.selected.(currHemi) = bi;
    pRFs.selected.(currHemi) = ci;
    commonCoords.(currHemi) = ai; 

    mismatchedDataVoxels = false(1,dataRuns.totalNrVoxels.(currHemi));
    mismatchedDataVoxels(~bi) = true;

    mismatchedPRFs = false(1,pRFs.totalNrVoxels.(currHemi));
    mismatchedPRFs(~ci) = true;

    dataRuns.mismatchedVoxels.(currHemi) = mismatchedDataVoxels;
    pRFs.mismatchedVoxels.(currHemi) = mismatchedPRFs;
        
end
    
out = struct();
out.pRFs = pRFs;
out.dataRuns = dataRuns;
out.coords = coords;
out.commonCoords = commonCoords;
