function [selectedVoxels, ve] = selectVoxelsFromPRFVarExpl(paramsROI, roiOrder, combROIs, thresh)

if isempty(thresh) || ~exist('thresh','var')
    thresh = 0.1;
end

% Get nr subjects and rois
[nrSubjects, nrROIs]  = size(paramsROI);

% Preallocate space
selectedVoxels = cell(nrSubjects,nrROIs);
ve = selectedVoxels;


% Loop over subjects
for sj = 1:nrSubjects
    
    % .. rois
    for idx = roiOrder
        curr_ve = [];
        % check if roi data exists
        if ~isempty(paramsROI{sj, roiOrder(idx)}) && ~ismember(roiOrder(idx),combROIs+1)
            
            % pRF var expl threshold to select voxels
            
            % Check if we want to combine this ROI (say VO1/VO2)
            if ismember(roiOrder(idx),combROIs) && ~isempty(paramsROI{sj, roiOrder(idx+1)})
                curr_ve = [paramsROI{sj, roiOrder(idx)}, paramsROI{sj, roiOrder(idx+1)}];
            end
            
            selectedVoxels{sj,roiOrder(idx)} = find(curr_ve>=thresh);
            if isempty(selectedVoxels{sj,roiOrder(idx)})
                selectedVoxels{sj,roiOrder(idx)} = NaN(size(curr_ve));
            end
            ve{sj,roiOrder(idx)} = curr_ve(selectedVoxels{sj,roiOrder(idx)});
            
        else
            curr_ve = paramsROI{sj, roiOrder(idx)};
            selectedVoxels{sj,roiOrder(idx)} = find(selectedVoxels{sj,roiOrder(idx)});
            ve{sj,roiOrder(idx)} = curr_ve(selectedVoxels{sj,roiOrder(idx)});
        end
        
    end
end