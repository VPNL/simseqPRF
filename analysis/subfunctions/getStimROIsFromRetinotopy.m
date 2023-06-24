function [] = getStimROIsFromRetinotopy(projectDir, subjnr, varargin)

% Function to get pRFs with center falling within stimulus corner and
% within 4 individual squares.
%
% Step 1. Run this function, save only large ROIS (saveLargeROIs=true, saveSmallROIs=false).
% Step 2. Check ROIs on indiv subject surface, and manually correct if needed.
% Step 3. Run this function again, now save only small ROIS (saveSmallROIs=true, saveLargeROIs=false).
% Step 4. Check ROIs again on indiv subject surface, and manually correct if needed.
%% 1. Set params
p = inputParser;
p.addRequired('projectDir',@ischar);
p.addRequired('subjnr', @isnumeric);
p.addParameter('hemis',{'lh','rh'}, @(x) any(validatestring(x{:},{'lh','rh'})));
p.addParameter('saveEntireCornerROIs', true, @islogical);
p.addParameter('saveIndivSquareROIs',false,@islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('saveFigs',true,@islogical);
p.addParameter('useStimFlag',true,@islogical);
p.addParameter('relativeSizes',{'small','large'}, ...
    @(x) any(validatestring(x{:},{'small','large'})));
p.parse(projectDir,subjnr,varargin{:});

% Rename variables
projectDir              = p.Results.projectDir;
subjnr                  = p.Results.subjnr;
hemis                   = p.Results.hemis; % Hemispheres to use
saveEntireCornerROIs    = p.Results.saveEntireCornerROIs; % Save large corner ROI, containing pRFs with centers overlapping 4 squares and space in between
saveIndivSquareROIs     = p.Results.saveIndivSquareROIs;  % Save individual square ROIs, only containing pRFs with centers overlapping single square
verbose                 = p.Results.verbose;  % If verbose is true, plot figures
saveFigs                = p.Results.saveFigs; % Save figure with pRFs
useStimFlag             = p.Results.useStimFlag; % if true, get stimulus locations and size from stimulus file, if false, use default simulation stim
relativeSizes           = p.Results.relativeSizes; %'small' or 'large'; % For pilot 3, varying square size

%% Set up paths
pths = getSubjectPaths(projectDir,subjnr);
cd(fullfile(pths.dataDirToon, pths.toon))

% Loop over hemispheres and sizes
for h = 1:length(hemis)
    for rs  = 1:length(relativeSizes)
        
        nrStimSquares = 4;          % 16; % choose 4 (pilot1 and 3) or 16 (pilt 2)
        
        if strcmp(relativeSizes{rs},'small')
            roiName_postFix = '_area2sq_eccen5';
        else
            roiName_postFix = '_area4sq_eccen5';
        end
        offset        = 0.41; % deg
        
        if useStimFlag,
            squareLoc = [];
            nrSquares = nrStimSquares;
        else
            if nrStimSquares==4
                squareLoc =  [0, 31, 32, 20, 21];  % 4 squares (each 1deg2)
            elseif nrStimSquares==16
                squareLoc = [0, 41:44, 30:33, 19:22, 8:11];  % 16 squares (each 1deg2)
            elseif nrStimSquares==9
                [71:2:75, 41:2:45, 11:2:15]; % Pilot 3: 9 squares (each 1.75 deg2)
            end
            nrSquares = length(squareLoc)-1;
        end
        
        % colors to plot pRFs for individual squares
        if isempty(nrSquares)
            colors = turbo(4+1);
        else
            colors = turbo(nrSquares+1);
        end
        
        if saveEntireCornerROIs || saveIndivSquareROIs
            if ~exist(fullfile('./3DAnatomy/ROIs/simseqROIs'),'dir')
                mkdir(fullfile('./3DAnatomy/ROIs/simseqROIs'));
            end
        end
        
        %% Get mrVista params and volume
        dt       = 'Averages'; % dataTYPE
        prfFit   = 'cssFit';   % prf model
        rois     =  pths.allROIs; % list of ROIs to load
        nRois    = length(rois);
        veThresh = 0.1; % 10% variance explained threshold
        
        % Get stimulus locations from file
        loc = getStimLocDeg(projectDir, subjnr, sessionNr, useStimFlag, squareLoc, hemis{h}, relativeSizes{rs}, offset);
        
        if strcmp(hemis{h}, 'rh')
            axeslim = [-12,2];
        elseif strcmp(hemis{h}, 'lh')
            axeslim = [-2,12];
        end
        
        if verbose
            makeprettyfigures;
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
        
        % Set T1 anatomy path and initialize hidden mrVista gray view
        setVAnatomyPath('./3DAnatomy/t1.nii.gz')
        hvol = initHiddenGray;
        hvol = viewSet(hvol, 'curdt', dt);
        
        %% Load results from pRF model
        
        rm_model = dir(fullfile('./Gray', dt, sprintf('*%s*fFit.mat',prfFit)));
        hvol = rmSelect(hvol,1,fullfile(rm_model.folder,rm_model.name));
        
        %% Load ROI coords and get pRF results
        
        if verbose
            fH1 = figure(1); clf; set(gcf, 'color', 'w', 'Position',[73,606,1994,415]);
            fH2 = figure(2); clf; set(gcf, 'color', 'w', 'Position',[73,606,1994,415]);
        end
        
        labels = []; for m=1:size(loc.xBounds,1), labels{m} = sprintf('square %d',m); end
        
        for r = 1:nRois
            
            % Get ROI name and path
            roiName = pths.preferredRoiName(hemis{h},rois{r});
            roiPath = fullfile('./3DAnatomy/ROIs/',sprintf('%s.mat',roiName));
            
            if exist(roiPath, 'file')
                
                % load ROI
                load(roiPath);
                
                % Get ROI and gray indices
                [~, indROI, indGray] = intersectCols(ROI.coords,hvol.coords);
                
                % Create same length ROI with NaNs, only put in data for indexed vertices
                roiData.name     = roiName;
                roiData.x0       = hvol.rm.retinotopyModels{1}.x0(indGray);
                roiData.y0       = hvol.rm.retinotopyModels{1}.y0(indGray);
                roiData.sigma    = hvol.rm.retinotopyModels{1}.sigma.major(indGray);
                roiData.exponent = hvol.rm.retinotopyModels{1}.exponent(indGray);
                if strcmp(prfFit, 'css')
                    roiData.effectiveSize     = roiData.sigma ./ sqrt(roiData.exponent);
                else
                    roiData.effectiveSize     = roiData.sigma;
                end
                radiusToPlot = 1; %2.5; % radius (in deg)
                %         roiData.sizeToPlot     = pi*(radiusToPlot*roiData.effectiveSize).^2;
                
                % Get variance explained and its ROI mask
                roiData.varexpl  = rmGet(hvol.rm.retinotopyModels{1}, 'varexplained');
                roiData.varexpl  = roiData.varexpl(indGray);
                roiData.veMask   = roiData.varexpl>veThresh;
                
                % Get polar angle and eccentricity
                [ang,ecc] = cart2pol(roiData.x0,roiData.y0);
                ang_deg = rad2deg(ang); clear ang;
                roiData.ecc = ecc;
                roiData.ang = ang_deg;
                
                clear square;
                if useStimFlag
                    x0_stim =  (roiData.x0>=minX) & (roiData.x0<=maxX);
                    y0_stim =  (roiData.y0>=minY) & (roiData.y0<=maxY);
                    totalFourStimROI_idx = sum([x0_stim;y0_stim;roiData.veMask],1)==3;
                    
                    for sq = 1:size(loc.xBounds,1)
                        square(:,sq) = (roiData.x0 >= loc.xBounds(sq,1)) & (roiData.x0 <= loc.xBounds(sq,2)) & (roiData.y0 >= loc.yBounds(sq,1)) & (roiData.y0 <= loc.yBounds(sq,2));
                    end
                else
                    if strcmp(hemis{h}, 'rh')
                        x0_stim =  (roiData.x0<=minX) & (roiData.x0>=maxX);
                        y0_stim =  (roiData.y0<=minY) & (roiData.y0>=maxY);
                        totalFourStimROI_idx = sum([x0_stim;y0_stim;roiData.veMask],1)==3;
                        
                        for sq = 1:size(loc.xBounds,1)
                            square(:,sq) = (roiData.x0 <= loc.xBounds(sq,1)) & (roiData.x0 >= loc.xBounds(sq,2)) & (roiData.y0 <= loc.yBounds(sq,1)) & (roiData.y0 >= loc.yBounds(sq,2));
                        end
                        
                    elseif strcmp(hemis{h}, 'lh')
                        x0_stim =  (roiData.x0>=minX) & (roiData.x0<=maxX);
                        y0_stim =  (roiData.y0>=minY) & (roiData.y0<=maxY);
                        totalFourStimROI_idx = sum([x0_stim;y0_stim;roiData.veMask],1)==3;
                        
                        for sq = 1:size(loc.xBounds,1)
                            square(:,sq) = (roiData.x0 >= loc.xBounds(sq,1)) & (roiData.x0 <= loc.xBounds(sq,2)) & (roiData.y0 >= loc.yBounds(sq,1)) & (roiData.y0 <= loc.yBounds(sq,2));
                        end
                    end
                end
                % Get indices within square and above variance expl threshold
                onlyInsideSquaresROI_idx = sum([square,roiData.veMask'],2)==2;
                
                % Copy and update ROI struct;
                ROIorig = ROI; clear ROI;
                
                % save large ROI
                ROI = ROIorig;
                ROI.coords = ROIorig.coords(:,indROI(totalFourStimROI_idx));
                ROI.created = datestr(now);
                ROI.modified = [];
                ROI.name = sprintf('%s_simseqstimcorner%d%s',roiName,nrSquares,roiName_postFix);
                ROI.roiSize = sprintf('%d voxels, %dmm^3',size(ROI.coords,2),size(ROI.coords,2));
                if saveEntireCornerROIs
                    save(fullfile('./3DAnatomy/ROIs/simseqROIs',sprintf('%s.mat',ROI.name)),'ROI');
                end
                clear ROI;
                
                % save semi-large ROI (only centers within squares)
                ROI = ROIorig;
                ROI.coords = ROIorig.coords(:,indROI(onlyInsideSquaresROI_idx));
                ROI.created = datestr(now);
                ROI.modified = [];
                ROI.name = sprintf('%s_simseqsquaresonly%d%s',roiName,nrSquares,roiName_postFix);
                ROI.roiSize = sprintf('%d voxels, %dmm^3',size(ROI.coords,2),size(ROI.coords,2));
                if saveEntireCornerROIs
                    save(fullfile('./3DAnatomy/ROIs/simseqROIs',sprintf('%s.mat',ROI.name)),'ROI');
                end
                clear ROI;
                
                % save small single square ROIs from simseqstimcorner ROI
                for ii = 1:size(square,2)
                    ROI = ROIorig;
                    ROI.coords = ROIorig.coords(:,indROI(square(:,ii)));
                    ROI.created = datestr(now);
                    ROI.modified = [];
                    ROI.name = sprintf('%s_simseqsquare%d_%d%s',roiName,ii,nrSquares,relativeSizes{rs});
                    ROI.roiSize = sprintf('%d voxels, %dmm^3',size(ROI.coords,2),size(ROI.coords,2));
                    if saveIndivSquareROIs
                        save(fullfile('./3DAnatomy/ROIs/simseqROIs',sprintf('%s.mat',ROI.name)),'ROI');
                    end
                    clear ROI;
                    
                    % Plot all and selected pRF XY centers
                    if verbose
                        
                        nrRows = 2;
                        nrCols = ceil(nRois/2);
                        figure(fH1);
                        subplot(nrRows,nrCols,r); hold all
                        set(gca,'Units','Normalized')
                        if ii == 1
                            plot([-25,25], [0,0],'k', [0,0],[-25,25], 'k'); hold all;
                        end
                        axis equal;  xlim(axeslim); ylim(axeslim)
                        %                 scatter(roiData.x0(square(:,ii)),roiData.y0(square(:,ii)),roiData.sizeToPlot(square(:,ii)), colors(ii,:))
                        for jj = find(square(:,ii))'
                            rectangle('Position',[roiData.x0(jj)-((radiusToPlot*roiData.effectiveSize(jj))/2),...
                                roiData.y0(jj)-((radiusToPlot*roiData.effectiveSize(jj))/2),...
                                radiusToPlot*roiData.effectiveSize(jj),radiusToPlot*roiData.effectiveSize(jj)],...
                                'Curvature',[1 1], ...
                                'EdgeColor',colors(ii,:));
                        end
                        
                        set(gca, 'FontSize', 15, 'TickDir', 'out')
                        if r>nrCols; xlabel('x position (deg)'); end
                        if (r==1 || r==nrCols+1) ylabel('y position (deg)'); end
                        title(strrep(roiName,'_', ' '))
                    end
                end
                %         if r==1
                %             l = findobj(gca);
                %             h = legend(l([(1:length(labels))+1]),labels,'location','none'); pos = get(h,'position');
                %             set(h, 'position',[0.9 0.45 0.1 0.1]); legend boxoff;
                %         end
                
                
                % Plot all and selected pRF XY centers
                if verbose
                    nrRows = 2;
                    nrCols = ceil(nRois/2);
                    figure(fH2);
                    subplot(nrRows,nrCols,r)
                    plot([-25,25], [0,0],'k', [0,0],[-25,25], 'k'); hold all;
                    for ii = find(roiData.veMask)
                        rectangle('Position',[roiData.x0(ii)-((radiusToPlot*roiData.effectiveSize(ii))/2),...
                            roiData.y0(ii)-((radiusToPlot*roiData.effectiveSize(ii))/2),...
                            radiusToPlot*roiData.effectiveSize(ii),radiusToPlot*roiData.effectiveSize(ii)], 'Curvature',[1 1], ...
                            'EdgeColor',[0.9290    0.6940    0.1250]);
                    end
                    %             scatter(roiData.x0(roiData.veMask),roiData.y0(roiData.veMask),roiData.sizeToPlot(roiData.veMask))
                    scatter(roiData.x0(totalFourStimROI_idx),roiData.y0(totalFourStimROI_idx),'rx')
                    scatter(roiData.x0(onlyInsideSquaresROI_idx),roiData.y0(onlyInsideSquaresROI_idx),'cx')
                    %% debug indices
                    %             [~, ~, subindices2] = intersectCols(ROI_stim.coords,hvol.coords);
                    %             scatter(hvol.rm.retinotopyModels{1}.x0(subindices2),hvol.rm.retinotopyModels{1}.y0(subindices2),'gx')
                    %         polarplot(roiData.ang(roiData.veMask),roiData.ecc(roiData.veMask),'o'); hold all;
                    %         polarplot(roiData.ang(stim_roi_idx), roiData.ecc(stim_roi_idx),'rx')
                    
                    xlim([-25,25]); ylim([-25,25])
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
            else
                warning('[%s]: Cannot find ROI: %s, skipping.. \n', mfilename, roiPath);
            end
        end
        
        if verbose && saveFigs
            if ~exist(fullfile(pths.figureDir, 'prfCenters'),'dir'); mkdir(fullfile(pths.figureDir, 'prfCenters')); end
            print(fH1,fullfile(pths.figureDir, 'prfCenters', sprintf('pilot3_%s_allpRFcenters_and_withinStim_%s_%d%s', hemis{h},prfFit,nrSquares,relativeSizes{rs})), '-depsc')
            print(fH2,fullfile(pths.figureDir, 'prfCenters', sprintf('pilot3_%s_pRFcentersWithinIndivSquares_%s_%d%s', hemis{h},prfFit,nrSquares,relativeSizes{rs})), '-depsc')
        end
    end
end
