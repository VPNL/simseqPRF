function fH = plotDataModelFitDoG_singleVoxel_StimulusConditionBlocks(projectDir,subjnrs,varargin)

% Parse inputs
p = inputParser;
p.addRequired('projectDir', @ischar); % '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal';
p.addRequired('subjnrs', @isnumeric); % Subject nrs are: [1,2,3,7,8,9,10,11,12,13]
p.addParameter('plotModelFlag',true, @islogical);
p.addParameter('roisToPlot',[],@(x) (iscell(x) || isnumeric(x)));
p.addParameter('roiType','stimcorner4_area4sq_eccen5',@ischar);
p.addParameter('hemi','both',@(x) any(validatestring(x,{'lh','rh','both'})));
p.addParameter('selectedDataVoxels',[],@isnumeric);
p.addParameter('spatialModel',{'onegaussianFit','cssFit','onegaussianFit','differenceofgaussianFit'}, @iscell);
p.addParameter('temporalModel',{'1ch-glm','1ch-glm','3ch-stLN','1ch-glm'}, @iscell);
p.addParameter('mdllbl',{'LSS','CSS','CST','DoG'}, @iscell);
p.addParameter('nrbslTRs',4,@isnumeric);
p.addParameter('nrows',5, @isnumeric);
p.addParameter('ncols',8, @isnumeric);
p.addParameter('sbplOrder',[4,2,8,6,3,1,7,5], @isnumeric);
p.addParameter('xl',[-5 15], @isnumeric);
p.addParameter('maxTrial',19, @isnumeric);
p.addParameter('saveFigs',false, @islogical);
p.addParameter('saveFigDir',[],@ischar);
p.addParameter('subDir',[],@ischar);
p.parse(projectDir,subjnrs,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff



% Set default plotting params
makeprettyfigures;

% Set matching Sim-Seq conditions
matchedCondsSmall = [1,5;2,6];
matchedCondsLarge = [3,7;4,8];
matchedConds      = [matchedCondsSmall;matchedCondsLarge];
stimlift          = 0.4; % Offset from x-axis for example stimulus sequence

% Define condition labels
subTtl = {'SEQ small 0.2s';'SEQ small 1s'; ...
    'SEQ large 0.2s';'SEQ large 1s'; ...
    'SIM small 0.2s';'SIM small 1s'; ...
    'SIM large 0.2s';'SIM large 1s'};

% Define colors
cmap = [getColormapPRFModels(0); ...
        getColormapPRFModels(2)];

if plotModelFlag
    colLinModel = cmap(1,:);
    colCSSModel = cmap(2,:);
    colCTSModel = cmap(3,:);
    colDoGModel = cmap(4,:);
end

for subjnr  = subjnrs 
   
    % Set subject paths
    pths       = getSubjectPaths(projectDir,subjnr);
    if ~exist(fullfile(pths.simseqResultsDir),'dir'), 
        mkdir(fullfile(pths.simseqResultsDir));
    end
    cd(fullfile(pths.simseqResultsDir))
    dataFolder = fullfile(pths.dataDirSimSeq, pths.subjID);

    if saveFigs
        saveFigDir = fullfile(pths.figureDir);
    end
    
    % Get ROIs depending on selected hemi(s)
    if strcmp(hemi,'both')
        rois       = pths.definedROIsBOTH;
    elseif strcmp(hemi,'lh')
        rois       = pths.definedROIsLH;
    elseif strcmp(hemi,'rh')
        rois       = pths.definedROIsRH;
    end
    
    if ~isempty(roisToPlot)
        if isnumeric(roisToPlot)
            rois = rois(roisToPlot);
        elseif iscell(roisToPlot)
            rois = rois(ismember(rois,roisToPlot));
        end
    end
    
    % Get stimulus file
    fname    = sprintf('stim_simseq_run2_v%d.mat', ...
        pths.expversionNr);
    stimfile = fullfile(pths.stimDir, pths.subjID, fname);
    
    % Get both small and large stimulus to plot with pRF location
    a       = load(stimfile);
    idx     = find(a.sequence(2,:)==max(a.sequence(1,:)));
    stimL   = ~squeeze(a.images(:,:,idx(1)));
    idx     = find(a.sequence(2,:)==6);
    stimS   = ~squeeze(a.images(:,:,idx(1)));
    clear a idx
    
    % define figures
    fH = figure(102); clf;  szScr = get(0,'ScreenSize');
    if plotModelFlag
        set(fH, 'Position', [0,0,szScr(3),szScr(4)]); hold on;
    else
        set(fH, 'Position', [0,0,szScr(3),szScr(4)/2]); hold on;
    end
    
    %% Loop over selected rois
    for roi = 1:length(rois) %
        clear p

        % Load predictions with data
        for pp = 1:length(spatialModel)
            
            modelName = mdllbl{pp};

            switch modelName
                case {'LSS','CSS','CST'}
                    
                    if ~exist('stimMS','var')
                        % Load stim time course in millisecond
                        dataFileToLoad = fullfile(dataFolder, sprintf('preprocData_%s_%s_%s_run12.mat', ...
                            rois{roi},roiType,'cssFit'));
                        d = load(dataFileToLoad);
                        stimMS = d.stimMS;
                        clear d;
                    end 
                    
                    predFileToLoad = fullfile(dataFolder,sprintf('%s_%s_%s_%s_cvfitResults.mat', ...
                        pths.subjID, rois{roi}, spatialModel{pp}, temporalModel{pp}));
                    p.(modelName) = load(predFileToLoad, 'T', 'allDataTrials',...
                        'params','R2_crossval_mean','data','B_crossval_mean');
                    
                    % Clear some memory
                    p.(modelName).params.stim.images = [];
                    p.(modelName).params.stim.images_unconvolved = [];
                    
                    % Rename prf params
                    p.(modelName).ve    = p.(modelName).data.params.varexpl.both;
                    p.(modelName).x0    = p.(modelName).data.params.x0.both;
                    p.(modelName).y0    = p.(modelName).data.params.y0.both;
                    p.(modelName).szeff = p.(modelName).data.params.effectiveSize.both;
                    
                    p.(modelName).exp      = p.(modelName).data.params.exponent.both;
                    p.(modelName).simseqR2 = p.(modelName).R2_crossval_mean;
                    
                    if strcmp(modelName, 'CST')
                        p.(modelName).exp_temporal = p.(modelName).data.params.exponent_temporal.both;
                        p.(modelName).beta_s = p.(modelName).B_crossval_mean(:,1);
                        p.(modelName).beta_t = p.(modelName).B_crossval_mean(:,2);
                    end
            
            
                case 'DoG'
                    
                    % Load data and model predictions (stored in one)
                    dataFileToLoad = fullfile(dataFolder, 'DogSimulation', sprintf('dogSimulation_%s_%s_%s.mat', ...
                                pths.subjID, rois{roi},roiType));
                    d = load(dataFileToLoad);
                    
                    p.(modelName).T.Condition      = d.TDoG.Condition;
                    p.(modelName).T.ModelTSCVMn    = d.TDoG.TS;
                    p.(modelName).T.ModelTSCVError = d.TDoG.TSerror;
                    p.(modelName).allDataTrials    = d.allDoGTrials;
                    p.(modelName).data.params      = d.params;
                    p.(modelName).simseqR2         = d.R2_crossval_mean;
                    p.(modelName).B_crossval_mean  = squeeze(mean(d.B_crossval,3));
                    
                    clear d;
                    p.(modelName).x0    = p.(modelName).data.params.analysis.spatial.x0;
                    p.(modelName).y0    = p.(modelName).data.params.analysis.spatial.y0;
                    p.(modelName).szeff = p.(modelName).data.params.analysis.spatial.effectiveSize;
                    p.(modelName).exp     = p.(modelName).data.params.analysis.spatial.exponent;
                    p.(modelName).sigmaCenter = p.(modelName).data.params.analysis.spatial.sigmaMajor;
                    p.(modelName).sigmaSurround = p.(modelName).data.params.analysis.spatial.sigmaSurround;
            end
        end % MODELS
        
        % Get selection of voxels to plot, based on R2
        if isempty(selectedDataVoxels)
            [~,idx_R2] = sort(p.CST.simseqR2,'descend');
            if length(idx_R2)<50, topVox = length(idx_R2); else topVox = 50; end
            selectedDataVoxels = idx_R2(1:topVox);
        end
        
        for pp = 1:length(spatialModel)
            p.(mdllbl{pp}).selectedPRFVoxels = selectedDataVoxels;
        end
        
        %% Let's plot
        for vv = selectedDataVoxels
            
            mnDataDC = {}; seDataDC = {};
            for cond = 1:size(matchedConds,1)
                
                % Select voxel data for condition
                currDataSEQ   = cell2mat(p.(mdllbl{1}).T.DataTSCV(p.(mdllbl{1}).T.Condition==matchedConds(cond,1)));
                currDataSEQSE = cell2mat(p.(mdllbl{1}).T.DataTSCVError(p.(mdllbl{1}).T.Condition==matchedConds(cond,1)));
                currDataSIM   = cell2mat(p.(mdllbl{1}).T.DataTSCV(p.(mdllbl{1}).T.Condition==matchedConds(cond,2)));
                currDataSIMSE = cell2mat(p.(mdllbl{1}).T.DataTSCVError(p.(mdllbl{1}).T.Condition==matchedConds(cond,2)));
                
                % baseline correct
                baselineSEQ = mean([currDataSEQ(1:nrbslTRs,vv)],'omitnan');
                baselineSIM = mean([currDataSIM(1:nrbslTRs,vv)],'omitnan');
                mnDataDC{matchedConds(cond,1)} = currDataSEQ(:,vv) - baselineSEQ;
                seDataDC{matchedConds(cond,1)} = currDataSEQSE(:,vv);
                mnDataDC{matchedConds(cond,2)} = currDataSIM(:,vv) - baselineSIM;
                seDataDC{matchedConds(cond,2)} = currDataSIMSE(:,vv);
                
                if plotModelFlag
                    % Model 1
                    currModelSEQ    = cell2mat(p.(mdllbl{1}).T.ModelTSCVMn(p.(mdllbl{1}).T.Condition==matchedConds(cond,1)));
                    currModelSEQSE  = cell2mat(p.(mdllbl{1}).T.ModelTSCVError(p.(mdllbl{1}).T.Condition==matchedConds(cond,1)));
                    currModelSIM    = cell2mat(p.(mdllbl{1}).T.ModelTSCVMn(p.(mdllbl{1}).T.Condition==matchedConds(cond,2)));
                    currModelSIMSE  = cell2mat(p.(mdllbl{1}).T.ModelTSCVError(p.(mdllbl{1}).T.Condition==matchedConds(cond,2)));
                    
                    if vv <= size(currModelSEQ,2)
                        mnModel1{matchedConds(cond,1)} = currModelSEQ(:,vv);
                        seModel1{matchedConds(cond,1)} = currModelSEQSE(:,vv);
                        
                        mnModel1{matchedConds(cond,2)} = currModelSIM(:,vv);
                        seModel1{matchedConds(cond,2)} = currModelSIMSE(:,vv);
                    else
                        mnModel1{matchedConds(cond,1)} = NaN(19,1);
                        seModel1{matchedConds(cond,1)} = NaN(19,1);
                        mnModel1{matchedConds(cond,2)} = NaN(19,1);
                        seModel1{matchedConds(cond,2)} = NaN(19,1);
                    end
                    
                    if length(mdllbl)>1
                        % Model 2
                        currModelSEQ    = cell2mat(p.(mdllbl{2}).T.ModelTSCVMn(p.(mdllbl{2}).T.Condition==matchedConds(cond,1)));
                        currModelSEQSE  = cell2mat(p.(mdllbl{2}).T.ModelTSCVError(p.(mdllbl{2}).T.Condition==matchedConds(cond,1)));
                        currModelSIM    = cell2mat(p.(mdllbl{2}).T.ModelTSCVMn(p.(mdllbl{2}).T.Condition==matchedConds(cond,2)));
                        currModelSIMSE  = cell2mat(p.(mdllbl{2}).T.ModelTSCVError(p.(mdllbl{2}).T.Condition==matchedConds(cond,2)));
                        
                        mnModel2{matchedConds(cond,1)} = currModelSEQ(:,vv);
                        seModel2{matchedConds(cond,1)} = currModelSEQSE(:,vv);
                        
                        mnModel2{matchedConds(cond,2)} = currModelSIM(:,vv);
                        seModel2{matchedConds(cond,2)} = currModelSIMSE(:,vv);
                    end
                    
                    if length(mdllbl)>2
                        % Model 3
                        currModelSEQ    = cell2mat(p.(mdllbl{3}).T.ModelTSCVMn(p.(mdllbl{3}).T.Condition==matchedConds(cond,1)));
                        currModelSEQSE  = cell2mat(p.(mdllbl{3}).T.ModelTSCVError(p.(mdllbl{3}).T.Condition==matchedConds(cond,1)));
                        currModelSIM    = cell2mat(p.(mdllbl{3}).T.ModelTSCVMn(p.(mdllbl{3}).T.Condition==matchedConds(cond,2)));
                        currModelSIMSE  = cell2mat(p.(mdllbl{3}).T.ModelTSCVError(p.(mdllbl{3}).T.Condition==matchedConds(cond,2)));
                        if vv <= size(currModelSEQ,2)
                            mnModel3{matchedConds(cond,1)} = currModelSEQ(:,vv);
                            seModel3{matchedConds(cond,1)} = currModelSEQSE(:,vv);
                            
                            mnModel3{matchedConds(cond,2)} = currModelSIM(:,vv);
                            seModel3{matchedConds(cond,2)} = currModelSIMSE(:,vv);
                        else
                            mnModel3{matchedConds(cond,1)} = NaN(19,1);
                            seModel3{matchedConds(cond,1)} = NaN(19,1);
                            mnModel3{matchedConds(cond,2)} = NaN(19,1);
                            seModel3{matchedConds(cond,2)} = NaN(19,1);
                            p.(mdllbl{3}).x0(vv) = NaN;
                            p.(mdllbl{3}).y0(vv) = NaN;
                            p.(mdllbl{3}).szeff(vv) = NaN;
                            p.(mdllbl{3}).ve(vv) = NaN;
                        end
                    end
                    
                    if length(mdllbl)>3
                        % Model 4
                        currModelSEQ    = cell2mat(p.(mdllbl{4}).T.ModelTSCVMn(p.(mdllbl{4}).T.Condition==matchedConds(cond,1)));
                        currModelSEQSE  = cell2mat(p.(mdllbl{4}).T.ModelTSCVError(p.(mdllbl{4}).T.Condition==matchedConds(cond,1)));
                        currModelSIM    = cell2mat(p.(mdllbl{4}).T.ModelTSCVMn(p.(mdllbl{4}).T.Condition==matchedConds(cond,2)));
                        currModelSIMSE  = cell2mat(p.(mdllbl{4}).T.ModelTSCVError(p.(mdllbl{4}).T.Condition==matchedConds(cond,2)));
                        if vv <= size(currModelSEQ,2)
                            mnModel4{matchedConds(cond,1)} = currModelSEQ(:,vv);
                            seModel4{matchedConds(cond,1)} = currModelSEQSE(:,vv);
                            
                            mnModel4{matchedConds(cond,2)} = currModelSIM(:,vv);
                            seModel4{matchedConds(cond,2)} = currModelSIMSE(:,vv);
                        else
                            mnModel4{matchedConds(cond,1)} = NaN(19,1);
                            seModel4{matchedConds(cond,1)} = NaN(19,1);
                            mnModel4{matchedConds(cond,2)} = NaN(19,1);
                            seModel4{matchedConds(cond,2)} = NaN(19,1);
                            p.(mdllbl{4}).x0(vv) = NaN;
                            p.(mdllbl{4}).y0(vv) = NaN;
                            p.(mdllbl{4}).szeff(vv) = NaN;
                            p.(mdllbl{4}).ve(vv) = NaN;
                        end
                    end
                end
            end
            
            %% Plot mean and SEM trials (data and model predictions)
            figure(fH); clf; hold on;
            numConditions = size(mnDataDC,2);
            % SEQ Small 200ms, SEQ Small 1s, SEQ Large 200ms, SEQ Large 1s, ...
            % SIM Small 200ms, SIM Small 1s, SIM Large 200ms, SIM Large 1s
            for cc = 1:numConditions
                
                % Get condition data
                mns = mnDataDC{cc};
                sems = seDataDC{cc};
                sbpl = sbplOrder(cc);
                stimTrial = stimMS{cc}';
                
                if exist('mnModel1','var')
                    mnsMod1 = mnModel1{cc};
                    semsMod1 = seModel1{cc};
                end
                if exist('mnModel2','var')
                    mnsMod2 = mnModel2{cc};
                    semsMod2 = seModel2{cc};
                end
                if exist('mnModel3','var')
                    mnsMod3 = mnModel3{cc};
                    semsMod3 = seModel3{cc};
                end
                if exist('mnModel4','var')
                    mnsMod4 = mnModel4{cc};
                    semsMod4 = seModel4{cc};
                end

                % Get stim trial, convert to binary
                stimTrial(stimTrial==1)=NaN;
                stimTrial(stimTrial==7)=NaN;
                uniqueStim = unique(stimTrial);
                for sq = 1:length(uniqueStim)
                    stimTrial(stimTrial==uniqueStim(sq)) = sq/50;
                end
                
                % Define stimulus time in seconds
                stimTime = [0:1:length(stimTrial)-1]./60;
                
                % Define x and y-axes
                tTrial = -4:(length(mns)-5);
                yl = [-2 3.5];
                if min(mns-sems) < yl(1)
                    yl(1) = min(mns-sems)-2;
                end
                if max(mns+sems) > yl(2)
                    yl(2) = 2+max(mns+sems);
                end
                
                % DATA + Linear Model
                ax = subplot(nrows,ncols,sbpl+ncols);  hold on;
                plot(tTrial(tTrial<=maxTrial),zeros(1,length(tTrial(tTrial<=maxTrial))),'k', 'lineWidth',0.5); hold on;
                plot(stimTime,[stimTrial+yl(1)+stimlift], 'color', 'k', 'lineWidth',2);
                s = shadedErrorBar(tTrial(tTrial<=maxTrial),mns(tTrial<=maxTrial)', sems(tTrial<=maxTrial)',...
                    'lineProps',{'Color','k'});hold on;
                s.mainLine.Visible = 'off';
                
                if plotModelFlag % Add LSS model fit
                    plot(tTrial(tTrial<=maxTrial),mnsMod1(tTrial<=maxTrial),':','LineWidth',3.5,'Color',colLinModel);
                else
                    xlabel('Time (s)');
                end
                
                title(sprintf('%s',  subTtl{cc}),'FontSize',12);
                box off;
                if any(isnan(yl)); yl = [-1 1]; end
                ylim([yl(1),yl(2)]); xlim(xl)
                if sbpl==1
                    if plotModelFlag
                        ylabel({sprintf('%s R2 %2.1f%%',mdllbl{1},100*p.LSS.simseqR2(vv)),'BOLD (% change)'});
                    else
                        ylabel('BOLD (% change)');
                    end
                else
                    set(gca,'YTick', [], 'YColor',[1 1 1])
                end
                
                if plotModelFlag
                    % DATA + CSS Model
                    if isfield(p, 'CSS')
                        subplot(nrows,ncols,sbpl+(2*ncols));
                        plot(tTrial(tTrial<=maxTrial),zeros(1,length(tTrial(tTrial<=maxTrial))),'k', 'lineWidth',0.5); hold on;
                        plot(stimTime,[stimTrial+yl(1)+stimlift], 'color', 'k', 'lineWidth',2);
                        s = shadedErrorBar(tTrial(tTrial<=maxTrial),mns(tTrial<=maxTrial)', sems(tTrial<=maxTrial)',...
                            'lineProps',{'Color','k'});hold on;
                        s.mainLine.Visible = 'off';
                        plot(tTrial(tTrial<=maxTrial),mnsMod2(tTrial<=maxTrial),'Color',colCSSModel, 'LineWidth',3.5,'LineStyle',':');
                        box off;
                        if any(isnan(yl)); yl = [-1 1]; end
                        ylim([yl(1),yl(2)]); xlim(xl)
                        if sbpl==1
                            ylabel({sprintf('%s R2 %2.1f%%',mdllbl{2},100*p.CSS.simseqR2(vv)), 'BOLD (% change)'});
                        else
                            set(gca,'YTick', [], 'YColor',[1 1 1])
                        end
                    end

                    % DATA + CST Model
                    if isfield(p, 'CST')
                        subplot(nrows,ncols,sbpl+(3*ncols));
                        plot(tTrial(tTrial<=maxTrial),zeros(1,length(tTrial(tTrial<=maxTrial))),'k', 'lineWidth',0.5); hold on;
                        plot(stimTime,[stimTrial+yl(1)+stimlift], 'color', 'k', 'lineWidth',2);
                        s = shadedErrorBar(tTrial(tTrial<=maxTrial),mns(tTrial<=maxTrial)', sems(tTrial<=maxTrial)',...
                            'lineProps',{'Color','k'});hold on;
                        s.mainLine.Visible = 'off';
                        plot(tTrial(tTrial<=maxTrial),mnsMod3(tTrial<=maxTrial),'Color',colCTSModel, 'LineWidth',3.5,'LineStyle',':');hold on;
                        box off;
                        if any(isnan(yl)); yl = [-1 1]; end
                        ylim([yl(1),yl(2)]); xlim(xl)

                        if ismember(sbpl, [1:ncols])
                            xlabel('Time (s)');
                        end
                        if sbpl==1
                            ylabel({sprintf('%s R2 %2.1f%%',mdllbl{3},100*p.CST.simseqR2(vv)), 'BOLD (% change)'});
                        else
                            set(gca,'YTick', [], 'YColor',[1 1 1])
                        end
                    end
                    
                    % DATA + DOG Model
                    if isfield(p, 'DoG')
                        subplot(nrows,ncols,sbpl+(4*ncols));
                        plot(tTrial(tTrial<=maxTrial),zeros(1,length(tTrial(tTrial<=maxTrial))),'k', 'lineWidth',0.5); hold on;
                        plot(stimTime,[stimTrial+yl(1)+stimlift], 'color', 'k', 'lineWidth',2);
                        s = shadedErrorBar(tTrial(tTrial<=maxTrial),mns(tTrial<=maxTrial)', sems(tTrial<=maxTrial)',...
                            'lineProps',{'Color','k'});hold on;
                        s.mainLine.Visible = 'off';
                        plot(tTrial(tTrial<=maxTrial),mnsMod4(tTrial<=maxTrial),'Color',colDoGModel, 'LineWidth',3.5,'LineStyle',':');hold on;
                        box off;
                        if any(isnan(yl)); yl = [-1 1]; end
                        ylim([yl(1),yl(2)]); xlim(xl)

                        if ismember(sbpl, [1:ncols])
                            xlabel('Time (s)');
                        end
                        if sbpl==1
                            ylabel({sprintf('%s R2 %2.1f%%',mdllbl{4},100*p.DoG.simseqR2(vv)), 'BOLD (% change)'});
                        else
                            set(gca,'YTick', [], 'YColor',[1 1 1])
                        end
                    end
                end
            end
            
            % Plot PRF + stim
            for fi = [1,2]
                if fi==1 % squareSize = 'Small'; 
                    stim = stimS;
                    loc = 1;
                else % squareSize = 'Large'; 
                    stim = stimL;
                    loc = 5;
                end
                ax1 = subplot(nrows,ncols,loc);
                thispRF.x0              = p.DoG.x0(vv);
                thispRF.y0              = p.DoG.y0(vv);
                thispRF.effectiveSize   = p.DoG.szeff(vv);
                thispRF.sigmaMajor      = p.DoG.sigmaCenter(vv);
                thispRF.sigmaSurround   = p.DoG.sigmaSurround(vv);
                thispRF.ve              = p.CST.ve(vv);
                thispRF.exp_spatial     = p.DoG.exp(vv);
                thispRF.beta            = p.DoG.B_crossval_mean(vv);
                ax1 = plotpRFlocOnImage(ax1,thispRF, stim, [], 'both');
                title(sprintf('pRF %d:[%2.1f,%2.1f], sz_c=%2.2f, sz_s=%2.2f, ve=%2.1f, beta=%2.1f', ...
                    vv, thispRF.x0,thispRF.y0,thispRF.sigmaMajor,...
                    thispRF.sigmaSurround,thispRF.ve,thispRF.beta), 'FontSize',9)
                xlabel(''); set(gca,'FontSize',9)
            end

            if plotModelFlag 
                sgtitle({sprintf('%s: Mean trial amplitude Data w/ Modelfit', rois{roi}), ''},...
                     'FontSize', 12,  'HorizontalAlignment','left')
                subplots = [(ncols+1):(ncols*nrows)];
            else
                sgtitle({sprintf('%s: Mean trial amplitude Data', rois{roi}), ''},...
                     'FontSize', 12,  'HorizontalAlignment','left')
                subplots = [ncols+(1:ncols)];
            end
            clear subYLim;
            
            for ll = 1:length(subplots)
                axTmp= subplot(nrows,ncols,subplots(ll));
                subYLim(ll,:) = axTmp.YLim;
            end
            minmaxrange = [min(subYLim(:,1)), max(subYLim(:,2))];
            sp = (min(subplots)-1) + find(subYLim(:,2)~=minmaxrange(2));
            if ~isempty(sp)
                for ll = sp'
                    axTmp = subplot(nrows,ncols,ll);
                    axTmp.YLim = minmaxrange;
                end
            end
            
            if saveFigs
               if ~exist('saveFigDir','var') || isempty(saveFigDir)
                    saveFigDir = fullfile(pths.figureDir);
                end
                subDir = 'supplFig7';
                subSaveFolder = fullfile(saveFigDir, subDir);
                if ~exist(subSaveFolder,'dir'); mkdir(subSaveFolder); end
                printName = sprintf('%s_singleVoxelTimeSeriesDoG_vox%d', ...
                    rois{roi},vv);
                print(gcf, fullfile(subSaveFolder,printName), '-dpng')
%                 print(gcf,'-depsc2','-painters','-r300','-loose',fullfile(subSaveFolder,printName));
            end
            
        end % VOXELS
    end % ROIS
end % SUBJECTS