function [selectedVoxels, rNC] = selectVoxelsFromNoiseCeiling(noiseceilingRun, amplTrial, roiOrder, combROIs, topX)

% Get nr subjects and rois
[nrSubjects, nrROIs]  = size(noiseceilingRun);

% Preallocate space
selectedVoxels = cell(nrSubjects,nrROIs);
rNC = selectedVoxels;

switch topX
    case 'top250vox'
        topN = 250;
        useTopX = true;
    case 'top50prct'
        prct = 0.5;
        useTopX = true;
end

% Loop over subjects
for sj = 1:nrSubjects
    
    % .. rois
    for idx = roiOrder
        
        % check if roi data exists
        if ~isempty(noiseceilingRun{sj, roiOrder(idx)}) && ~ismember(roiOrder(idx),combROIs+1)
            
            % Use 250 voxels with highest highest noise ceiling (i.e., splithalf reliability) 
            if useTopX && ~isempty(noiseceilingRun{sj, roiOrder(idx)})
                rNC{sj,roiOrder(idx)} = noiseceilingRun{sj, roiOrder(idx)};
                
                % Check if we want to combine this ROI (say VO1/VO2)
                if ismember(roiOrder(idx),combROIs) && ~isempty(noiseceilingRun{sj, roiOrder(idx+1)})
                    rNC{sj,roiOrder(idx)} = [noiseceilingRun{sj, roiOrder(idx)}, noiseceilingRun{sj, roiOrder(idx+1)}];
                end
                
                [~,maxNCRun] = sort(rNC{sj,roiOrder(idx)},'descend');
                if exist('prct','var')
                    topN = ceil(prct*length(maxNCRun));
                end
                if length(maxNCRun)<topN
                    maxVox = length(maxNCRun);
                else
                    maxVox = topN;
                end
                
                selectedVoxels{sj,roiOrder(idx)} = unique([maxNCRun(1:maxVox)]);
                
            else
                if ismember(roiOrder(idx),combROIs)
                    rNC{sj,roiOrder(idx)} = [noiseceilingRun{sj, roiOrder(idx)}, noiseceilingRun{sj, roiOrder(idx+1)}];
                    selectedVoxels{sj,roiOrder(idx)} = 1:(size(amplTrial{sj, 1, roiOrder(idx)},2) ...
                                            +size(amplTrial{sj, 1, roiOrder(idx+1)},2));
                else
                    selectedVoxels{sj,roiOrder(idx)} = 1:size(amplTrial{sj, 1, roiOrder(idx)},2);
                    rNC{sj,roiOrder(idx)} = noiseceilingRun{sj, roiOrder(idx)};
                end
            end
        end
    end
end