function fH = plotDataModelFitSTRet_singleVoxel_StimulusConditionBlocks(projectDir,subjnrs,varargin)

% Parse inputs
p = inputParser;
p.addRequired('projectDir', @ischar); % '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal';
p.addRequired('subjnrs', @isnumeric); % Subject nrs are: [1,2,3,7,8,9,10,11,12,13]
p.addParameter('plotModelFlag',true, @islogical);
p.addParameter('roisToPlot',[],@(x) (iscell(x) || isnumeric(x)));
p.addParameter('roiType','stimcorner4_area4sq_eccen5_stRet_CSTopt_DNST_matchingVoxels',@ischar);
p.addParameter('hemi','both',@(x) any(validatestring(x,{'lh','rh','both'})));
p.addParameter('selectedDataVoxels',[],@isnumeric)
p.addParameter('spatialModel',{'onegaussianFit','onegaussianFit','onegaussianFit'}, @iscell);
p.addParameter('temporalModel',{'3ch-stLN','3ch-stLN','1ch-dcts'}, @iscell);
p.addParameter('mdllbl',{'CSTfix','CSTopt','DNST'}, @iscell);
p.addParameter('nrbslTRs',4,@isnumeric);                                    % figure params: nr of baseline TRs, prior to block onset
p.addParameter('nrows',4, @isnumeric);                                      % figure params: nr of rows
p.addParameter('ncols',8, @isnumeric);                                      % figure params: nr of columns
p.addParameter('sbplOrder',[4,2,8,6,3,1,7,5], @isnumeric);                  % figure params: horizontal plotting order of conditions
p.addParameter('xl',19, @isnumeric);                                        % figure params: x-axis length (nr time points
p.addParameter('maxTrial',[-5 15], @isnumeric);                             % figure params: time axis of trial (seconds)
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

% Matched Sim-Seq conditions
matchedCondsSmall = [1,5;2,6];
matchedCondsLarge = [3,7;4,8];
matchedConds      = [matchedCondsSmall;matchedCondsLarge];

subTtl = {'SEQ small 0.2s';'SEQ small 1s'; ...
    'SEQ large 0.2s';'SEQ large 1s'; ...
    'SIM small 0.2s';'SIM small 1s'; ...
    'SIM large 0.2s';'SIM large 1s'};

xl = [-5 15];
maxTrial = 19;
stimlift = 0.4; % Offset from x-axis for example stimulus sequence

% Define colors
% tmp = getColormapPRFModels(0);

if plotModelFlag
    cmap = getColormapPRFModels(3);
end

for subjnr  = subjnrs

    % Set subject paths
    pths       = getSubjectPaths(projectDir,subjnr);
    if ~exist(fullfile(pths.simseqResultsDir),'dir'), 
        mkdir(fullfile(pths.simseqResultsDir));
    end
    cd(fullfile(pths.simseqResultsDir))
    dataFolderBase = fullfile(pths.dataDirSimSeq, pths.subjID);
    
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
    
    % Set up figure
    fH = figure(102); clf;  szScr = get(0,'ScreenSize');
    set(fH, 'Position', [0,0,szScr(3),szScr(4)]); hold on;
    
    %% Loop over selected rois
    for roi = 1:length(rois) %
        clear p
        % Load predictions with data
        for pp = 1:length(mdllbl)
            modelName = mdllbl{pp};
            dataFolder1 = sprintf('%s_preprocData_stRet_matchingVoxels', pths.subjID);
            
            switch modelName
                case {'CSTopt','DNST'}
                    if strcmp(modelName, 'CSTopt')
                        dataFolder2 = sprintf('preprocData_stRet_matchingVoxels_%s_opt',temporalModel{pp});
                    end
                    
                    % Load preprocessed time courses
                    dd = dir(fullfile(dataFolderBase, dataFolder1, ...
                        dataFolder2, 'varyStimDur*', sprintf('preprocData_%s_%s_%s_run12.mat', ...
                        rois{roi},roiType, spatialModel{pp})));
                    dataFileToLoad = fullfile(dd(1).folder,dd(1).name);
                    predFileToLoad = fullfile(dataFolderBase,...
                        sprintf('%s_%s_cvfits', pths.subjID,modelName),...
                        sprintf('%s_%s_%s_cvfits.mat', pths.subjID,rois{roi},modelName));
                   
                    p.(modelName) = load(predFileToLoad);
                    if strcmp('CSTopt',modelName)
                        p.(modelName).exp_temporal = p.(modelName).data.params.exponent_temporal.both;
                        p.(modelName).beta_s = p.(modelName).B_crossval_mean(:,1);
                        p.(modelName).beta_t = p.(modelName).B_crossval_mean(:,2);
                        p.(modelName).tau = p.(modelName).data.params.temporal.param.tau_s;
                    elseif strcmp('DNST',modelName)
                        p.(modelName).n = p.(modelName).params.analysis.temporal.param.n;
                        p.(modelName).semisat = p.(modelName).params.analysis.temporal.param.sigma;
                        p.(modelName).tau1    = p.(modelName).params.analysis.temporal.param.tau1;
                        p.(modelName).tau2    = p.(modelName).params.analysis.temporal.param.tau2;
                        p.(modelName).beta    = p.(modelName).B_crossval_mean(:,1);
                    end
                    
                case 'CSTfix'
                    dataFolder2 = sprintf('preprocData_stRet_matchingVoxels_%s_opt',temporalModel{pp});
                    dd = dir(fullfile(dataFolderBase, dataFolder1, dataFolder2, ...
                        'varyStimDur*', sprintf('preprocData_%s_%s_%s_run12.mat', ...
                        rois{roi},roiType,spatialModel{pp})));
                    dataFileToLoad = fullfile(dd(1).folder,dd(1).name);
                    predFileToLoad = fullfile(dataFolderBase,...
                        sprintf('%s_%s_cvfits', pths.subjID,modelName),...
                        sprintf('%s_%s_%s_cvfitResults.mat', pths.subjID,rois{roi},modelName));
                    
                    p.(modelName) = load(predFileToLoad);
                    
                    p.(modelName).exp_temporal = p.(modelName).data.params.exponent_temporal.both;
                    p.(modelName).beta_s = p.(modelName).B_crossval_mean(:,1);
                    p.(modelName).beta_t = p.(modelName).B_crossval_mean(:,2);
                    p.(modelName).tau = NaN(size(p.(modelName).beta_t));
            end
                  
            d = load(dataFileToLoad);
            stimMS = d.stimMS;
            clear d;

            % Clear some memory
            p.(modelName).params.stim.images = [];
            p.(modelName).params.stim.images_unconvolved = [];
            p.(modelName).ve       = p.(modelName).data.params.varexpl.both;
            p.(modelName).x0       = p.(modelName).data.params.x0.both;
            p.(modelName).y0       = p.(modelName).data.params.y0.both;
            p.(modelName).szeff    = p.(modelName).data.params.effectiveSize.both;
            p.(modelName).exp      = p.(modelName).data.params.exponent.both;
            p.(modelName).simseqR2 = p.(modelName).R2_crossval_mean;
            
        end % MODELS
        
        % Get selection of voxels to plot, based on R2
        if isempty(selectedDataVoxels)
            tmp = p.CSTopt.simseqR2;
            tmp(isnan(tmp))=0;
            [~,idx_R2] = sort(tmp,'descend');
            if length(idx_R2)<200, topVox = length(idx_R2); else topVox = 200; end
            selectedDataVoxels = idx_R2;%(1:topVox);
%             selectedDataVoxels = selectedDataVoxels(selectedDataVoxels<length(p.DN_ST.simseqR2));
        end
        
        for pp = 1:length(spatialModel)
            p.(mdllbl{pp}).selectedPRFVoxels = selectedDataVoxels;
        end
        
        %% Let's plot
        for vv = selectedDataVoxels
            
            mnDataDC = {}; seDataDC = {};
            for cond = 1:size(matchedConds,1)

                if plotModelFlag
                    for pp = 1:length(mdllbl)
                        % Select voxel data for condition
                        currDataSEQ   = cell2mat(p.(mdllbl{pp}).T.DataTSCV(matchedConds(cond,1)));
                        currDataSEQSE = cell2mat(p.(mdllbl{pp}).T.DataTSCVError(matchedConds(cond,1)));
                        currDataSIM   = cell2mat(p.(mdllbl{pp}).T.DataTSCV(matchedConds(cond,2)));
                        currDataSIMSE = cell2mat(p.(mdllbl{pp}).T.DataTSCVError(matchedConds(cond,2)));
                        
                        % baseline correct
                        baselineSEQ = mean([currDataSEQ(1:nrbslTRs,vv)],'omitnan');
                        baselineSIM = mean([currDataSIM(1:nrbslTRs,vv)],'omitnan');
                        mnDataDC{pp,matchedConds(cond,1)} = currDataSEQ(:,vv) - baselineSEQ;
                        seDataDC{pp,matchedConds(cond,1)} = currDataSEQSE(:,vv);
                        mnDataDC{pp,matchedConds(cond,2)} = currDataSIM(:,vv) - baselineSIM;
                        seDataDC{pp,matchedConds(cond,2)} = currDataSIMSE(:,vv);
                        
                        % Model 1
                        currModelSEQ    = cell2mat(p.(mdllbl{pp}).T.ModelTSCVMn(matchedConds(cond,1)));
                        currModelSEQSE  = cell2mat(p.(mdllbl{pp}).T.ModelTSCVError(matchedConds(cond,1)));
                        currModelSIM    = cell2mat(p.(mdllbl{pp}).T.ModelTSCVMn(matchedConds(cond,2)));
                        currModelSIMSE  = cell2mat(p.(mdllbl{pp}).T.ModelTSCVError(matchedConds(cond,2)));
                       
                        
                        if vv <= size(currModelSEQ,2)
                            mnModel{pp,matchedConds(cond,1)} = currModelSEQ(:,vv);
                            seModel{pp,matchedConds(cond,1)} = currModelSEQSE(:,vv);
                            
                            mnModel{pp,matchedConds(cond,2)} = currModelSIM(:,vv);
                            seModel{pp,matchedConds(cond,2)} = currModelSIMSE(:,vv);
                        else
                            mnModel{pp,matchedConds(cond,1)} = NaN(19,1);
                            seModel{pp,matchedConds(cond,1)} = NaN(19,1);
                            mnModel{pp,matchedConds(cond,2)} = NaN(19,1);
                            seModel{pp,matchedConds(cond,2)} = NaN(19,1);
                            p.(mdllbl{pp}).x0(vv) = NaN;
                            p.(mdllbl{pp}).y0(vv) = NaN;
                            p.(mdllbl{pp}).szeff(vv) = NaN;
                            p.(mdllbl{pp}).ve(vv) = NaN;
                        end
                        
                    end
                    
                    
                end
            end
            
            %% Plot mean and SEM trials (data and model predictions)
            figure(fH); clf; hold on;
            numConditions = size(mnDataDC,2);
            % SEQ Small 200ms, SEQ Small 1s, SEQ Large 200ms, SEQ Large 1s, ...
            % SIM Small 200ms, SIM Small 1s, SIM Large 200ms, SIM Large 1s
            for pp = 1:length(mdllbl)
                for cc = 1:numConditions
                    
                    % Get condition data
                    sbpl = sbplOrder(cc);
                    stimTrial = stimMS{cc}';
                    
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
                    tTrial = -4:(length(mnDataDC{pp,cc})-5);
                    yl = [-2 3.5];
                    if min(mnDataDC{pp,cc}-seDataDC{pp,cc}) < yl(1)
                        yl(1) = min(mnDataDC{pp,cc}-seDataDC{pp,cc})-2;
                    end
                    if max(mnDataDC{pp,cc}-seDataDC{pp,cc}) > yl(2)
                        yl(2) = 2+max(mnDataDC{pp,cc}+seDataDC{pp,cc});
                    end
                    
                    % DATA + Linear Model
                    ax = subplot(nrows,ncols,sbpl+(pp*ncols));  hold on;
                    plot(tTrial(tTrial<=maxTrial),zeros(1,length(tTrial(tTrial<=maxTrial))),'k', 'lineWidth',0.5); hold on;
                    plot(stimTime,[stimTrial+yl(1)+stimlift], 'color', 'k', 'lineWidth',2);
                    s = shadedErrorBar(tTrial(tTrial<=maxTrial),mnDataDC{pp,cc}(tTrial<=maxTrial)', seDataDC{pp,cc}(tTrial<=maxTrial)',...
                        'lineProps',{'Color','k'});hold on;
                    s.mainLine.Visible = 'off';
                    
                    if plotModelFlag % Add LSS model fit
                        plot(tTrial(tTrial<=maxTrial),mnModel{pp,cc}(tTrial<=maxTrial),':','LineWidth',3.5,'Color',cmap(pp,:));
                    else
                        xlabel('Time (s)');
                    end
                    
                    title(sprintf('%s',  subTtl{cc}),'FontSize',12);
                    box off;
                    if any(isnan(yl)); yl = [-1 1]; end
                    ylim([yl(1),yl(2)]); xlim(xl)
                    if sbpl==1
                        if plotModelFlag
                            ylabel({sprintf('%s R2 %2.1f%%',mdllbl{pp},100*p.(mdllbl{pp}).simseqR2(vv)),'BOLD (% change)'});
                        else
                            ylabel('BOLD (% change)');
                        end
                    else
                        set(gca,'YTick', [], 'YColor',[1 1 1])
                    end

                end
                
                
                stim = stimL+stimS;
                loc = 1 + ((pp-1)*3);

                ax1 = subplot(nrows,ncols,loc);
                thispRF.x0 = p.(mdllbl{pp}).x0(vv);
                thispRF.y0 = p.(mdllbl{pp}).y0(vv);
                thispRF.effectiveSize = p.(mdllbl{pp}).szeff(vv);
                thispRF.varexpl = p.(mdllbl{pp}).ve(vv);
                ax1 = plotpRFlocOnImage(ax1,thispRF, stim, [], 'both');
                if ismember(mdllbl{pp},{'CSTopt','CSTfix'})
                    thispRF.exp    = p.(mdllbl{pp}).exp_temporal(vv);
                    thispRF.tau    = p.(mdllbl{pp}).tau(vv);
                    thispRF.beta_s = p.(mdllbl{pp}).beta_s(vv);
                    thispRF.beta_t = p.(mdllbl{pp}).beta_t(vv);
                    title(sprintf('pRF %d:[%2.1f,%2.1f], sz=%2.2f, ve=%2.1f, exp=%1.2f, tau=%3.1f, beta=[%2.1f,%2.1f]', ...
                    vv, thispRF.x0,thispRF.y0,thispRF.effectiveSize,...
                    thispRF.varexpl, thispRF.exp, thispRF.tau*10, thispRF.beta_s,thispRF.beta_t), 'FontSize',9)
                else
                    thispRF.n       = p.(mdllbl{pp}).n(vv);
                    thispRF.semisat   = p.(mdllbl{pp}).semisat(vv);
                    thispRF.tau1    = p.(mdllbl{pp}).tau1(vv);
                    thispRF.tau2    = p.(mdllbl{pp}).tau2(vv);
                    thispRF.beta    = p.(mdllbl{pp}).B_crossval_mean(vv);
                    title(sprintf('pRF %d:[%2.1f,%2.1f], sz=%2.2f, ve=%2.1f,n=%2.2f,sigma=%2.2f, tau1=%3.1f,tau2=%3.1f,beta=[%2.1f]', ...
                    vv, thispRF.x0,thispRF.y0,thispRF.effectiveSize,...
                    thispRF.varexpl, thispRF.semisat,thispRF.semisat,thispRF.tau1,thispRF.tau2, thispRF.beta), 'FontSize',9)
                end
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
                subFigDir = 'supplfig8';
                subSaveFolder = fullfile(saveFigDir, subFigDir);
                if ~exist(subSaveFolder,'dir'); mkdir(subSaveFolder); end
                printName = sprintf('%s_singleVoxelTimeSeries_vox%d', ...
                    rois{roi},vv);
                print(gcf, fullfile(subSaveFolder,printName), '-dpng')
%                 print(gcf,'-depsc2','-painters','-r300','-loose',fullfile(subSaveFolder,printName));
            end
            
        end % VOXELS
    end % ROIS
end % SUBJECTS