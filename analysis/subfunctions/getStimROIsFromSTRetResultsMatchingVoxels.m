function [] = getStimROIsFromSTRetResultsMatchingVoxels(projectDir, subjnr, varargin)

% Function to get pRFs with center falling within stimulus corner and
% within 4 individual squares.
%
% Run this function, save only large ROIS (saveLargeROIs=true, saveSmallROIs=false).
% Check ROIs on indiv subject surface, and manually correct if needed.

%% 1. Set params
p = inputParser;
p.addRequired('projectDir',@ischar);
p.addRequired('subjnr', @isnumeric);
p.addParameter('temporalModels', {'3ch-stLN','1ch-dcts'}, @iscell);
p.addParameter('hemis',{'lh','rh'}, @(x) any(validatestring(x{:},{'lh','rh'})));
p.addParameter('saveEntireCornerROIs', true, @islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('saveFigs',true,@islogical);
p.addParameter('useStimFlag',true,@islogical);
p.parse(projectDir,subjnr,varargin{:});

% Rename variables
projectDir              = p.Results.projectDir;
subjnr                  = p.Results.subjnr;
temporalModels          = p.Results.temporalModels;
hemis                   = p.Results.hemis; % Hemispheres to use
saveEntireCornerROIs    = p.Results.saveEntireCornerROIs; % Save large corner ROI, containing pRFs with centers overlapping 4 squares and space in between
verbose                 = p.Results.verbose;  % If verbose is true, plot figures
saveFigs                = p.Results.saveFigs; % Save figure with pRFs

% ROI params
nrSquares       = 4;          % 16; % choose 4 (pilot1 and 3) or 16 (pilt 2)
roiType         = 'simseqstimcorner4_area4sq_eccen5';
roiName_postFix = ['_' roiType '_stRet_CSTopt_DNST_matchingVoxels'];
offset          = 0.41; % deg
colors          = turbo(nrSquares+1); % colors to plot pRFs for individual squares

if verbose
    makeprettyfigures;
end

%% Set up paths
sessionNr = getSessionNrMainExp(subjnr);
pths = getSubjectPaths(projectDir,subjnr,sessionNr);

cd(fullfile(pths.stRetDataDir, pths.stRet))

if saveEntireCornerROIs
    if ~exist(fullfile('./3DAnatomy/ROIs/simseqROIs'),'dir')
        mkdir(fullfile('./3DAnatomy/ROIs/simseqROIs'));
    end
end

for t = 1:length(temporalModels)
    % Rename temporalModel to match data table variables
    tmodel_full{t}     = strrep(temporalModels{t},'-','_');
    tmp             = strsplit(temporalModels{t},'-');
    tmodel{t}  = tmp{2};
end

%% Get stRet params and cortical surface volume
rm_stRet    = load(fullfile(pths.stRetResultsDir,'df_cv0_defaultHRF.mat'));
rm_stRet.DT.y0 = -1.* rm_stRet.DT.y0; % Known bug, stRet y-coords are flipped.
DT_simseq   = cell(1,length(temporalModels));

rois        =  pths.allROIs(1:find(ismember(pths.allROIs,{'IPS1'}))); % list of ROIs to load
nRois       = length(rois);
veThresh    = 0.1; % 10% variance explained threshold
useStimFlag = true;
radiusToPlot = 1; %2.5; % radius (in deg)

% Loop over hemispheres and sizes
for h = 1:length(hemis)
    
    % Get stimulus locations from file
    loc = getStimLocDeg(projectDir, subjnr, sessionNr, useStimFlag, nrSquares, hemis{h}, 'big', offset);
    
    if strcmp(hemis{h}, 'rh')
        axeslim = [-12,2];
    elseif strcmp(hemis{h}, 'lh')
        axeslim = [-2,12];
    end
    
    if useStimFlag && strcmp(hemis{h}, 'rh')
        loc.yBounds = -1.*loc.yBounds;
    end
    
    minX = min(loc.xBounds(:));
    maxX = max(loc.xBounds(:));
    minY = min(loc.yBounds(:));
    maxY = max(loc.yBounds(:));
    
    if ~useStimFlag && strcmp(hemis{h}, 'rh')
        minX = -minX;
        minY = -minY;
        maxX = -maxX;
        maxY = -maxY;
        
        loc.xBounds = -1.*loc.xBounds;
        loc.yBounds = -1.*loc.yBounds;
    end
    
    %     %% Load results from toonotopy pRF model
    %     hvol = viewSet(hvol, 'curdt', dt);
    %     prfFit = 'onegaussianFit';
    %     rm_model = dir(fullfile('./Gray', dt, sprintf('*%s*fFit.mat',prfFit)));
    %     hvol = rmSelect(hvol,1,fullfile(rm_model.folder,rm_model.name));
    
    %% Load ROI coords and get pRF results
    if verbose
        fH1 = figure(1); clf; set(gcf, 'color', 'w', 'Position',[73,606,1994,415]);
        fH2 = figure(2); clf; set(gcf, 'color', 'w', 'Position',[73,606,1994,415]);
    end
    
    labels = []; for m=1:size(loc.xBounds,1), labels{m} = sprintf('square %d',m); end
    
    for r = 1:length(rois)
        
        clear hvol coordsIndex_stRet_loc coordsMaskstRet coordsIndex_simseq_abs coordsMaskSimseq;
        
        %  Go to stRet data session
        cd(fullfile(pths.stRetDataDir, pths.stRet))
        
        % Set T1 anatomy path and initialize hidden mrVista gray view
        setVAnatomyPath('./3DAnatomy/t1.nii.gz')
        hvol = initHiddenGray;
        
        % Get Toon ROI name and path
        prfRoiName = pths.preferredRoiName(hemis{h},rois{r});
        roiPathToon = fullfile('./3DAnatomy/ROIs/',sprintf('%s.mat',prfRoiName));
        roiPathSimSeq = fullfile('./3DAnatomy/ROIs/simseqROIs',sprintf('%s_%s.mat',prfRoiName, roiType));
        
        if exist(roiPathToon, 'file')
            
            roiName = [hemis{h} '_' rois{r}];
            
            for t = 1:length(temporalModels)
                idx{t} = ismember(rm_stRet.DT.subjID,pths.subjID) & ...
                    ismember(rm_stRet.DT.roiName,roiName) & ...
                    ismember(rm_stRet.DT.model,'st') & ...
                    ismember(rm_stRet.DT.tmodel, tmodel_full{t});
                
                % Get pRF params
                prfParamsStRetExp{t} = rm_stRet.DT(idx{t},:);
            end
            
            if length(temporalModels)>1
                % check coordinates:
                if isequal(rm_stRet.DT.coords(idx{1}),rm_stRet.DT.coords(idx{2}))
                    stRetcoords = rm_stRet.DT.coords(idx{1});
                else
                    error('[%s]: Coordinates don''t match between temporal models!',mfilename)
                end
            end
            
            if ~isempty(stRetcoords)
                 
                %% Within stRet mrVista session
                % Get coords idx for entire toon ROI vs entire hemi
                hvol = loadROI(hvol,fullfile(roiPathToon),[],[],1,0);
                [coordsIdx_stRet_rel, coords3D_stRet_rel] = roiIndices(hvol, hvol.ROIs(1).coords);
                
                if ~isempty(setdiff(coordsIdx_stRet_rel,stRetcoords))
                    [outlier01,outlier01_i] = setdiff(coordsIdx_stRet_rel,stRetcoords);
                    coordsIdx_stRet_rel(outlier01_i) = [];
                end
                if ~isempty(setdiff(stRetcoords,coordsIdx_stRet_rel))
                    [outlier02,outlier02_i] = setdiff(stRetcoords,coordsIdx_stRet_rel);
                     stRetcoords(outlier02_i) = [];
                end
                assert(isequal(stRetcoords,coordsIdx_stRet_rel)); % table coords are local coords
                
                % Get preserved/absolute coordinates for entire ROI
                % ROIs(1) = entire ROI, e.g., "V1_toon"
                % ROIs(2) = simseq stim ROI, e.g., "V1_toon_simseqstimcorner4_area4sq_eccen5"
                [coordsIdx_stRet_abs, coords3D_stRet_abs] = roiIndices(hvol, hvol.ROIs(1).coords,1);
                
                if exist('outlier01','var') && (length(coordsIdx_stRet_abs(~isnan(coordsIdx_stRet_abs)))~=length(coordsIdx_stRet_rel))
                    tmp0 = coordsIdx_stRet_abs(~isnan(coordsIdx_stRet_abs));
                    numIdx = find(~isnan(coordsIdx_stRet_abs));
                    outlier_idx = ismember(tmp0,outlier01);
                    coordsIdx_stRet_abs(numIdx(outlier_idx))=NaN;
                end
                
                % If equal, we have correspondance between relative and
                % preserved (absolute) coordinates:
                % where: coords3D_stRet_rel <==> coords3D_stRet_abs
                assert(isequal(coords3D_stRet_abs, hvol.ROIs(1).coords));
                
                % Step 1: Get preserved simseq stim corner ROI coords using stRet mrVista Session
                hvol = loadROI(hvol,fullfile(roiPathSimSeq),[],[],1,0);
                [coordsIdx_simseq_abs, coords3D_simseq_abs] = roiIndices(hvol, hvol.ROIs(2).coords,1);
                
                % Step 2: Check if table coordinates overlap with table
                coords3D_simseq_abs_noNan  = coords3D_simseq_abs(:,~isnan(coordsIdx_simseq_abs));
                coordsIdx_simseq_abs_noNan = coordsIdx_simseq_abs(~isnan(coordsIdx_simseq_abs));
                coords3D_stRet_abs_NoNan   = coords3D_stRet_abs(:,~isnan(coordsIdx_stRet_abs));
                coordsIdx_stRet_abs_NoNan  = coordsIdx_stRet_abs(~isnan(coordsIdx_stRet_abs));
                
                if ~isempty(coords3D_simseq_abs_noNan)
                    % MAIN STEP: get common coordinates stRet prfs with simseq corner
                    [commonCoords_stimcornerSimSeq_stRet,stimcornerIdx_simseq,stimcornerIdx_stRet] = ...
                        intersectCols(coords3D_simseq_abs_noNan,coords3D_stRet_abs_NoNan);
                    
                    % Put coords back in ROI array (no NaNs) -- 3x1626 vertices:
                    stRetCoords3D_within_simseqstimcorner_NoNan = NaN(size(coords3D_stRet_abs_NoNan));
                    stRetCoords3D_within_simseqstimcorner_NoNan(:,stimcornerIdx_stRet)= coords3D_stRet_abs_NoNan(:,stimcornerIdx_stRet);
                    assert(isequal(coords3D_stRet_abs_NoNan(:,stimcornerIdx_stRet), commonCoords_stimcornerSimSeq_stRet));
                    
                    % Also make a boolean -- 1626x1
                    stRetCoords3D_within_simseqstimcorner_bool = false(size(coords3D_stRet_abs_NoNan,2),1);
                    stRetCoords3D_within_simseqstimcorner_bool(stimcornerIdx_stRet)= true;
                    
                    % Put coords back in full array -- 3x2366 vertices
                    stRetCoords3D_within_simseqstimcorner = NaN(size(coords3D_stRet_abs));
                    stRetCoords3D_within_simseqstimcorner(:,~isnan(coordsIdx_stRet_abs)) = stRetCoords3D_within_simseqstimcorner_NoNan;
                    
                    % Get the actual coordinate numbers of the stRet table,
                    % and check if they correspond to simseq coord Indices
                    stRetCoords_label = stRetcoords(stRetCoords3D_within_simseqstimcorner_bool);
                    
                    coordsIdx_simseq_abs_nonNan2 = coordsIdx_simseq_abs_noNan(stimcornerIdx_simseq)';
                    
                    % check if there is any outliers for stRetCoords and SimSeq
                    % Coords
                    if ~isempty(setdiff(coordsIdx_simseq_abs_nonNan2,stRetCoords_label)) || ...
                            ~isempty(setdiff(stRetCoords_label,coordsIdx_simseq_abs_nonNan2))

                        stRetCoords_label2 = intersect(stRetCoords_label,coordsIdx_simseq_abs_nonNan2);
                        stRetCoords3D_within_simseqstimcorner_bool2 = stRetCoords3D_within_simseqstimcorner_bool;
                        [outlier1, outlier1_i] = setdiff(stRetCoords_label,stRetCoords_label2);
                        for ol = 1:length(outlier1)
                            xx = find(stRetCoords3D_within_simseqstimcorner_bool2);
                            if outlier1_i(ol)<=length(xx)
                                stRetCoords3D_within_simseqstimcorner_bool2(xx(outlier1_i(ol)))=false;
                                stRetCoords_label(outlier1_i(ol)) = [];
                                commonCoords_stimcornerSimSeq_stRet(:,outlier1_i(ol)) = [];
                            end
                        end
%                         threeDcoordOutlier = sum([stRetCoords3D_within_simseqstimcorner_bool2,stRetCoords3D_within_simseqstimcorner_bool],2)==1;

                        stRetCoords3D_within_simseqstimcorner_bool = stRetCoords3D_within_simseqstimcorner_bool2;
                    end
                                        
                    % Step 3: If overlapping, subselect prfParamsStRetExp.coords that fall within simseq stim corner ROI.
                    % so that we only use those preserved ROI coords that overlap with simseq stimulus corner ROI.
                    %                 [ai,a] = ismember(stRetCoordsIdx_within_simseqstimcorner,coordsIdx_simseq_abs_noNan);
                    %                 assert(~isempty(sum(ai)))
                    
                    % Step 4: check which table coords correspond to these coords:
                    %                 [bi,b] = ismember(rm_stRet.DT.coords(idx{1}),coordsIdx_stRet_abs_NoNan(stimcornerIdx_stRet));
                    %                 [ci,c] = ismember(rm_stRet.DT.coords(idx{2}),coordsIdx_stRet_abs_NoNan(stimcornerIdx_stRet));
                    
                    %                 [b,ci] = find(prfParamsStRetExp{1}.coords == stRetCoordsIdx_within_simseqstimcorner);
                    %                 assert(isequal(prfParamsStRetExp{1}.coords(b),stRetCoordsIdx_within_simseqstimcorner(ci)'));
                    %                 [b,ci] = find(prfParamsStRetExp{2}.coords == stRetCoordsIdx_within_simseqstimcorner);
                    %                 assert(isequal(prfParamsStRetExp{2}.coords(b),stRetCoordsIdx_within_simseqstimcorner(ci)'));
                    
                    %% Get final coordinates in 3D and indices
                    finalCoords3D = commonCoords_stimcornerSimSeq_stRet;
                    %                 finalCoordIdx{1} = find(bi); %coordsIdx_stRet_abs_NoNan(stimcornerIdx_stRet);
                    %                 finalCoordIdx{2} = find(ci);
                    %                 isequal(finalCoordIdx{1},finalCoordIdx{2});
                    
                    
                    % for stRet Exp. get bool indices
                    
                    selectedPRFs0 = [];
                    for t = 1:length(temporalModels)
                        x0       = prfParamsStRetExp{t}.x0(stRetCoords3D_within_simseqstimcorner_bool);
                        y0       = prfParamsStRetExp{t}.y0(stRetCoords3D_within_simseqstimcorner_bool);
                        varexpl  = prfParamsStRetExp{t}.cv_varexp(stRetCoords3D_within_simseqstimcorner_bool);
                        veMask   = varexpl>veThresh;
                        
                        x0_stim =  (x0>=minX) & (x0<=maxX);
                        y0_stim =  (y0>=minY) & (y0<=maxY);
                        totalFourStimROI_idx = sum([x0_stim,y0_stim,veMask],2)==3;
                        selectedPRFs0(:,t) = totalFourStimROI_idx;
                    end
                    
                    selectedPRFs = (sum(selectedPRFs0,2)==2);
                    
                    if ~isempty(selectedPRFs)
                        
                        % Copy and update ROI struct;
                        ROIorig = hvol.ROIs(2); clear hvol;
                        
                        % save ROI
                        ROI = ROIorig;
                        ROI.coords = finalCoords3D(:,selectedPRFs);
                        ROI.created = datestr(now);
                        ROI.modified = [];
                        ROI.comments = 'Toonotopy ROI border for vox within SimSeq stim corner, with stRet pRF params (defaultHRF); ROI coords stRet session matching';
                        ROI.name = sprintf('%s%s',prfRoiName,roiName_postFix);
                        ROI.roiSize = sprintf('%d voxels, %dmm^3',size(ROI.coords,2),size(ROI.coords,2));
                        save(fullfile('./3DAnatomy/ROIs/simseqROIs',sprintf('%s.mat',ROI.name)),'ROI');
                        
                        clear ROI
                        %% store pRF params
                        for t = 1:length(temporalModels)
                            
                            clear roiData
                            %                     % get surface coords temporal model:
                            %                     x = (rm_stRet.DT.coords(idx{t})); % 1 by 1626 vertices
                            %                     assert(isequal(x,coordsIdx_stRet_abs_NoNan'))
                            
                            % get surface coords matching across models, with
                            x2 = find(stRetCoords3D_within_simseqstimcorner_bool);
                            stRetCoords3D_within_simseqstimcorner_bool_matching = false(size(stRetCoords3D_within_simseqstimcorner_bool)); % 1 by 1626 vertices
                            stRetCoords3D_within_simseqstimcorner_bool_matching(x2(selectedPRFs)) = true;
                            
                            % check if subset of voxels overlaps in x0 values.
                            %                     isequal(ismember(prfParamsStRetExp{t}.x0(stRetCoords3D_within_simseqstimcorner_bool_matching),...
                            %                         prfParamsStRetExp{t}.x0(stRetCoords3D_within_simseqstimcorner_bool))
                            
                            % Create same length ROI with NaNs, only put in data for indexed vertices
                            roiData.name     = roiName;
                            roiData.x0       = prfParamsStRetExp{t}.x0(stRetCoords3D_within_simseqstimcorner_bool_matching);
                            roiData.y0       = prfParamsStRetExp{t}.y0(stRetCoords3D_within_simseqstimcorner_bool_matching);
                            roiData.sigma    = prfParamsStRetExp{t}.sigma(stRetCoords3D_within_simseqstimcorner_bool_matching);
                            roiData.temporal = prfParamsStRetExp{t}.temporal(stRetCoords3D_within_simseqstimcorner_bool_matching);
                            
                            isequal(prfParamsStRetExp{t}.coords(stRetCoords3D_within_simseqstimcorner_bool_matching), stRetCoords_label(selectedPRFs));
                            isequal(prfParamsStRetExp{t}.coords(stRetCoords3D_within_simseqstimcorner_bool_matching), coordsIdx_simseq_abs_noNan(selectedPRFs)')
                            roiData.coords   = prfParamsStRetExp{t}.coords(stRetCoords3D_within_simseqstimcorner_bool_matching);
                            
                            %                     roiData.effectiveSize     = roiData.sigma; % do not create effectiveSize yet
                            
                            % Get variance explained and its ROI mask
                            roiData.varexpl  = prfParamsStRetExp{t}.cv_varexp(selectedPRFs);
                            roiData.veMask   = roiData.varexpl>veThresh;
                            
                            % Get polar angle and eccentricity
                            [ang,ecc] = cart2pol(roiData.x0,roiData.y0);
                            ang_deg = rad2deg(ang); clear ang;
                            roiData.ecc = ecc;
                            roiData.ang = ang_deg;
                            
                            
                            prfParamsStRetExp0 = prfParamsStRetExp{t}(stRetCoords3D_within_simseqstimcorner_bool_matching,:);
                            prfParamsStRetExp1 = [prfParamsStRetExp0,table(finalCoords3D(:,selectedPRFs)','VariableNames',{'3DCoords'})];
                            
                            DT_simseq{t} = [DT_simseq{t}; prfParamsStRetExp1];
                            
                            clear prfParamsStRetExp0 prfParamsStRetExp1
                            
                            
                            % Plot all and selected pRF XY centers
                            if verbose
                                nrRows = 2;
                                nrCols = ceil(nRois/2);
                                figure(fH2);
                                subplot(nrRows,nrCols,r)
                                plot([-25,25], [0,0],'k', [0,0],[-25,25], 'k'); hold all;
                                for ii = find(roiData.veMask')
                                    rectangle('Position',[roiData.x0(ii)-((radiusToPlot*roiData.sigma(ii))/2),...
                                        roiData.y0(ii)-((radiusToPlot*roiData.sigma(ii))/2),...
                                        radiusToPlot*roiData.sigma(ii),radiusToPlot*roiData.sigma(ii)], 'Curvature',[1 1], ...
                                        'EdgeColor',colors(t,:));
                                end
                                %                         scatter(roiData.x0,roiData.y0,'rx')
                                %% debug indices
                                %             [~, ~, subindices2] = intersectCols(ROI_stim.coords,hvol.coords);
                                %             scatter(hvol.rm.retinotopyModels{1}.x0(subindices2),hvol.rm.retinotopyModels{1}.y0(subindices2),'gx')
                                %         polarplot(roiData.ang(roiData.veMask),roiData.ecc(roiData.veMask),'o'); hold all;
                                %         polarplot(roiData.ang(stim_roi_idx), roiData.ecc(stim_roi_idx),'rx')
                                
                                xlim([-20,20]); ylim([-20,20])
                                axis square; set(gca, 'FontSize', 15, 'TickDir', 'out')
                                if r==1 || r ==nrCols+1
                                    ylabel('y position (deg)');
                                end
                                if r>nrCols
                                    xlabel('x position (deg)');
                                end
                                title(strrep(roiName,'_', ' '))
                                %             if r==nRois
                                %                 l = findobj(gca);
                                %                 h = legend(l([4,3,2]),{sprintf('all pRF centers (ve>%d%%)',veThresh*100), 'within corner', 'within square'}, 'Location', 'None');
                                %                 set(h, 'position',[0.85 0.3 0.1 0.1]); legend boxoff;legend boxoff;
                                %             end
                            end
                            
                        end
                        
                    else
                        warning('[%s]: No overlap with ROI: %s, skipping.. \n', mfilename, roiPathToon);
                    end
                else
                    warning('[%s]: No overlap with ROI: %s, skipping.. \n', mfilename, roiPathToon);
                end
                
            else
                warning('[%s]: No stRet ROI data: %s, skipping.. \n', mfilename, roiPathToon);
            end
        else
            warning('[%s]: Cannot find ROI: %s, skipping.. \n', mfilename, roiPathToon);
        end
    end
    
    if verbose && saveFigs
        if ~exist(fullfile(pths.figureDir, 'pRFs'),'dir'); mkdir(fullfile(pths.figureDir, 'pRFs')); end
        print(fH1,fullfile(pths.figureDir, 'pRFs', sprintf('pRFCenters_stRet_CSToptDNST_matchingVoxels_%s_allpRFcenters_and_withinStim_%d%s', hemis{h},nrSquares,'big')), '-depsc')
        print(fH2,fullfile(pths.figureDir, 'pRFs', sprintf('pRFCenters_stRet_CSToptDNST_matchingVoxels_%s_pRFcentersWithinIndivSquares_%d%s', hemis{h},nrSquares,'big')), '-depsc')
    end
end

DT_simseq_names = temporalModels;
save(fullfile(pths.stRetResultsDir,sprintf('%s_stRet_DT_CSTopt_DNST_matchingVoxels_simseqSession.mat',pths.subjID)),'DT_simseq','DT_simseq_names')

