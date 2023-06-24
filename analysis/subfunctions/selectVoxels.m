function [selectedVoxels, rNC, ve] = selectVoxels(noiseceilingRun, paramsROI, nc_thresh, ve_thresh, roiOrder, combROIs)

% Get nr subjects and rois
[nrSubjects, nrROIs]  = size(noiseceilingRun);

% Preallocate space
selectedVoxels = cell(nrSubjects,nrROIs);
rNC = selectedVoxels;
ve  = selectedVoxels;

% Loop over subjects
for sj = 1:nrSubjects
    
    % .. rois
    for idx = 1:length(roiOrder)
        
        % check if roi data exists and if we want to combine this ROI (say VO1/VO2)
        if ~isempty(noiseceilingRun{sj, roiOrder(idx)}) && ismember(roiOrder(idx),combROIs)
            
            % Threshold based on splithalf reliability
            if isempty(paramsROI{sj, roiOrder(idx+1)}) && isempty(noiseceilingRun{sj, roiOrder(idx+1)})
                ve{sj,roiOrder(idx)}  = [paramsROI{sj, roiOrder(idx)}.varexpl.both];
                rNC{sj,roiOrder(idx)} = [noiseceilingRun{sj, roiOrder(idx)}];
            else
                ve{sj,roiOrder(idx)}  = [paramsROI{sj, roiOrder(idx)}.varexpl.both, paramsROI{sj,roiOrder(idx+1)}.varexpl.both];
                rNC{sj,roiOrder(idx)} = [noiseceilingRun{sj, roiOrder(idx)}, noiseceilingRun{sj, roiOrder(idx+1)}];
            end
        else 
            if ~isempty(noiseceilingRun{sj, roiOrder(idx)}) && ~ismember(roiOrder(idx),combROIs+1)
                rNC{sj,roiOrder(idx)} = noiseceilingRun{sj, roiOrder(idx)};
                ve{sj,roiOrder(idx)}  = paramsROI{sj, roiOrder(idx)}.varexpl.both;
            else
                rNC{sj,roiOrder(idx)} = [];
                ve{sj,roiOrder(idx)}  = [];
            end
%             selectedVoxels{sj,roiOrder(idx)} = 1:(size(amplTrial{sj, 1, roiOrder(idx)},2) ...
%                 +size(amplTrial{sj, 1, roiOrder(idx+1)},2));
                
%                     selectedVoxels{sj,roiOrder(idx)} = 1:size(amplTrial{sj, 1, roiOrder(idx)},2);
%                     rNC{sj,roiOrder(idx)} = noiseceilingRun{sj, roiOrder(idx)};
           
        end
        if ~isempty(rNC{sj,roiOrder(idx)}) && ~isempty(ve{sj,roiOrder(idx)})
            selectedVoxels{sj,roiOrder(idx)} = ((rNC{sj,roiOrder(idx)}>=nc_thresh) & (ve{sj,roiOrder(idx)}>=ve_thresh));
        else
            selectedVoxels{sj,roiOrder(idx)} = [];
        end
    end
end