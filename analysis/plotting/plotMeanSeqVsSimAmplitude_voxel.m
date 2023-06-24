function fH = plotMeanSeqVsSimAmplitude_voxel(ds, fLMM, lmmResults, varargin)
% Function to plot average voxel SEQ block amplitude vs SIM block
% amplitude, colored by pRF size

% Parse inputs
p = inputParser;
p.addRequired('ds');
p.addRequired('fLMM', @iscell);
p.addRequired('lmmResults', @isstruct);
p.addParameter('LMMlbl',[],@ischar);
p.addParameter('subjnrs',[1:3,7:13],@isnumeric);
p.addParameter('roisToPlot',[], @isnumeric);
p.addParameter('roiType','stimcorner4_area4sq_eccen5',@ischar);
p.addParameter('conditionNamesSimSeq',{'Small-0.2s','Small-1s','Big-0.2s','Big-1s'}, @iscell);
p.addParameter('conditionOrderSimSeq',[1:4],@isnumeric);
p.addParameter('ve_thresh',0.2, @isnumeric); % threshold in percentage noise ceiling from simseq exp
p.addParameter('nc_thresh',0.1, @isnumeric); % threshold in percentage of pRF var expl. from toonotopy
p.addParameter('plotAllSubjectsTogether',false, @islogical);
p.addParameter('plotModelAmpl',false, @islogical);
p.addParameter('whichModelToPlot','LSS', @ischar);
p.addParameter('plotDataFlag',true, @islogical);
p.addParameter('plotFitFlag',true, @islogical);
p.addParameter('plotModelAmplSubj',false, @islogical);
p.addParameter('saveFigs',true, @islogical);
p.addParameter('saveFigDir',[],@ischar);
p.addParameter('subDir',[],@ischar);
p.addParameter('mrkrSz',24, @isnumeric);
p.addParameter('edgesSz',[0:2:12], @isnumeric); % Color dots by pRF size (binned)
p.addParameter('AlphaLevelMarker',1, @isnumeric);
p.parse(ds,fLMM,lmmResults,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff

%%
allRoisToPlot = unique(ds.ROI,'stable');
newROIOrder = [1,2,3,4,5,8,9,6,7];
allRoisToPlot = allRoisToPlot(newROIOrder);
if isempty(roisToPlot)
    roisToPlot = allRoisToPlot;
elseif ~isempty(roisToPlot)
    if isnumeric(roisToPlot)
        tmp = allRoisToPlot(roisToPlot);
        roisToPlot = tmp;
    elseif iscell(roisToPlot)
        tmp = allRoisToPlot(ismember(allRoisToPlot,roisToPlot));
        roisToPlot = tmp;
    end
end

if isempty(saveFigDir)
    saveFigDir = fullfile(simseqRootPath, 'results','group');
end

if plotModelAmplSubj
    cmapSubjects = gray(10);
end

% Plotting params
cmapSz   = parula(length(edgesSz));
cmapSz   = cmapSz(1:end-1,:);
plotOrder = [2,1,4,3];

%% Choose all vs single subject
if plotAllSubjectsTogether
    
    axRange = [-1,7];
    if ~isfield(lmmResults,'fixedIntercepts') || ~isfield(lmmResults,'fixedSlopes') || ...
            ~isfield(lmmResults,'fixedIntercepts_CI') || ~isfield(lmmResults,'fixedSlopes_CI')
        error('[%s]: Group Slope and intercept w/ CI need to be an input variable', mfilename)
    end
    
    if length(conditionOrderSimSeq)==1
        fH = figure; set(gcf,'Position',[60, 438, 1920, 474]); %[ 1  1 1920 950]);
        ncols = length(roisToPlot);
    end
    for idx = 1:length(roisToPlot)
        if length(conditionOrderSimSeq)>1
            fH = figure(50); clf; set(gcf,'Position',[60, 438, 1920, 474]); %[ 1  1 1920 950]);
            ncols = length(conditionOrderSimSeq);
        end
        
        % Get PRF CST size
        szCSTToPlot = ds.pRFCSTsize(...
            ds.ROI==nominal(roisToPlot(idx)) & ...
            ds.Condition==nominal(1));
        
        for c = conditionOrderSimSeq
            if plotModelAmpl
                fn_seq = ['MeanSeqAmpModel' whichModelToPlot];
                fn_sim = ['MeanSimAmpModel' whichModelToPlot];
                % Get SEQ and SIM amplitude model prediction
                xToPlot = ds.(fn_seq)(...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
                yToPlot = ds.(fn_sim)(...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
            else
                % Get SEQ and SIM amplitude data
                xToPlot = ds.MeanSeqAmp(...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
                yToPlot = ds.MeanSimAmp(...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
            end
            if ~isempty(xToPlot)
                if length(conditionOrderSimSeq)==1
                    subplot(1,ncols,idx); hold on;
                    xlabel('SEQ (% signal)','FontSize',8);
                    if idx==1, ylabel('SIM (% signal)','FontSize',8); end
                else
                    subplot(1,ncols,plotOrder(c)); hold on;
                    ylabel('SIM (% signal)','FontSize',8);
                    xlabel('SEQ (% signal)','FontSize',8);
                end
                title(sprintf('%s: %s',string(roisToPlot(idx)),  conditionNamesSimSeq{c}),'FontSize',10);

                plot([-5 21],[0 0],'k','lineWidth',0.25)
                plot([0 0],[-5 21],'k','lineWidth',0.25)
                plot([-5 21],[-5 21],'k:','lineWidth',0.25)
                set(gca, 'FontSize',10);
                box off; axis square
                set(gca,'xlim', axRange, 'ylim',axRange);
                
                if plotDataFlag
                    [~,vIdx] = sort(szCSTToPlot,'ascend');
                    bb = discretize(szCSTToPlot(vIdx),edgesSz);
                    scatter(xToPlot(vIdx),yToPlot(vIdx),mrkrSz, szCSTToPlot(vIdx),'filled',...
                        'MarkerEdgeColor',[1 1 1],'MarkerEdgeAlpha',AlphaLevelMarker,...
                        'MarkerFaceAlpha',AlphaLevelMarker); hold on;
                    cb = colorbar; cb.Position = cb.Position + 0.01;
                    colormap(cmapSz);
                    set(gca, 'CLim', [edgesSz(1),edgesSz(end)]);
                    cb.Label.String = 'effective pRF size (deg)';
                    cb.TickDirection = 'out'; cb.Box = 'off';
                    cb.Ticks = edgesSz;
                end
                
                if plotFitFlag
                    xPoints = [min(xToPlot),max(xToPlot)];
                    yfixed = lmmResults.fixedIntercepts{idx}(c) + xPoints.*lmmResults.fixedSlopes{idx}(c);
                    yfixedCI_Lower = lmmResults.fixedIntercepts_CI{idx}(c,1) + xPoints.*lmmResults.fixedSlopes_CI{idx}(c,1);
                    yfixedCI_Upper = lmmResults.fixedIntercepts_CI{idx}(c,2) + xPoints.*lmmResults.fixedSlopes_CI{idx}(c,2);
               
                    
                    plot([min(xToPlot),max(xToPlot)],yfixed,'k','lineWidth',2)
                    plot([min(xToPlot),max(xToPlot)],yfixedCI_Lower,'k:','lineWidth',2)
                    plot([min(xToPlot),max(xToPlot)],yfixedCI_Upper,'k:','lineWidth',2)
                end
                
                if plotModelAmplSubj
                    for ss = 1:10
                        voxIdx = (tmpT.Subject==ss);
                        if ~isempty(voxIdx)
                            predYSubj{ss} = {xToPlot(voxIdx), yfixed(voxIdx)};
                            plot(xToPlot(voxIdx),ypred(voxIdx),':','color',cmapSubjects(ss,:),'lineWidth',2); hold on;
                        end
                    end
                end
            end
        end
        
        if (ve_thresh > 0) || (nc_thresh > 0)
            if plotModelAmpl
                sgtitle(sprintf('All Subjects: Seq vs Sim predicted amplitudes by %s (voxels with pRF ve>%1.1f & nc>%1.1f)', ...
                    whichModelToPlot, ve_thresh, nc_thresh))
            else
                sgtitle(sprintf('All Subjects: Seq vs Sim amplitudes (voxels with pRF ve>%1.1f & nc>%1.1f)', ...
                    ve_thresh, nc_thresh))
            end
        else
            if plotModelAmpl
                sgtitle(sprintf('All Subjects: Seq vs Sim predicted amplitudes by %s (all voxels)', ...
                    whichModelToPlot))
            else
                sgtitle('All Subjects: Seq vs Sim amplitudes (all voxels)')
            end
        end
        
        % Save figure
        if saveFigs
            if plotModelAmpl
                fName = sprintf('Summary%sAllSubjects_SEQvsSIM_%sModelAmpl_XVal1_%s',...
                    replace(string(roisToPlot(idx)),'/',''), whichModelToPlot, roiType);
            else
                fName = sprintf('Summary%sAllSubjects_SEQvsSIMAmpl_XVal1_%s',...
                    replace(string(roisToPlot(idx)),'/',''), roiType);
            end
            if plotFitFlag
                fName = sprintf('%s_LMMfit_%s',fName, LMMlbl);
                if plotDataFlag
                    fName = [fName '_wData'];
                else
                    fName = [fName '_noData'];
                end
            else
                fName = sprintf('%s_dataOnly',fName);
            end
            
            thisSaveFigDir = fullfile(saveFigDir);
            if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
            saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%             print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
        end
    end
    
    %% Plot subjects separately
else
    axRange  = [[-1,5];[-1,5];[-1,7]; ...
        [-1,7];[-1,7];[-1,7]; ...
        [-1,4];[-1,4]; [-1,4]]; ...
        
    % Loop over subjects
    for sj = 1:length(subjnrs)
        % Loop over ROIs
        for idx = 1:length(roisToPlot)
            fH = figure(50); clf; set(gcf,'Position',[60, 438, 1920, 474]);
            
            for c = conditionOrderSimSeq
                % Get PRF CST size and exp
                szCSTToPlot = ds.pRFCSTsize(ds.Subject==subjnrs(sj) & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
                
                % Get SEQ and SIM amplitude data
                xToPlot = ds.MeanSeqAmp(ds.Subject==subjnrs(sj) & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
                yToPlot = ds.MeanSimAmp(ds.Subject==subjnrs(sj) & ...
                    ds.ROI==nominal(roisToPlot(idx)) & ...
                    ds.Condition==nominal(c));
                
                % Create temporary data table to get single subject fit
                tmpT = table();
                tmpT.MeanSeqAmp = xToPlot;
                tmpT.MeanSimAmp = yToPlot;
                tmpT.Subject    = ones(size(xToPlot)).*(subjnrs(sj));
                tmpT.ROI        = repmat(roisToPlot(idx),size(xToPlot));
                tmpT.Condition  = ones(size(xToPlot)).*c;
                
                tmpT           = table2dataset(tmpT);
                tmpT.Subject   = tmpT.Subject;
                tmpT.Condition = nominal(tmpT.Condition);
                tmpT.ROI       = nominal(tmpT.ROI);
                
                if ~isempty(xToPlot)
                    [ypred, ypredCI] = predict(fLMM{(allRoisToPlot==roisToPlot(idx))}, tmpT);
                    
                    subplot(1,4,plotOrder(c)); hold on;
                    title(sprintf('%s: %s',string(unique(tmpT.ROI)),  conditionNamesSimSeq{c}),'FontSize',10);
                    ylabel('SIM (% signal)','FontSize',10);
                    xlabel('SEQ (% signal)','FontSize',8);
                    
                    plot([-5 21],[0 0],'k','lineWidth',0.25)
                    plot([0 0],[-5 21],'k','lineWidth',0.25)
                    plot([-5 21],[-5 21],'k:','lineWidth',0.25)
                    set(gca, 'FontSize',10);
                    box off; axis square
                    set(gca,'xlim', axRange(idx,:), 'ylim',axRange(idx,:));
                    
                    if plotDataFlag
                        [~,vIdx] = sort(szCSTToPlot,'ascend');
                        scatter(xToPlot(vIdx),yToPlot(vIdx), mrkrSz, szCSTToPlot(vIdx), 'filled',...
                            'MarkerEdgeColor',[1 1 1],'MarkerEdgeAlpha',AlphaLevelMarker,...
                            'MarkerFaceAlpha',AlphaLevelMarker); hold on;
                        
                        cb = colorbar;
                        colormap(cmapSz);
                        set(gca, 'CLim', [edgesSz(1),edgesSz(end)]);
                        cb.Label.String = 'effective pRF size (deg)';
                        cb.TickDirection = 'out'; cb.Box = 'off';
                        cb.Ticks = edgesSz;
                        
                    end
                    
                    if plotFitFlag
                        plot(xToPlot,ypred,':','color','k','lineWidth',1.5); hold on;
                    end
                end
            end
            
            if (ve_thresh > 0) || (nc_thresh > 0)
                if plotModelAmpl
                    sgtitle(sprintf('Subject S%d: Seq vs Sim predicted amplitudes by %s (voxels with pRF ve>%1.1f & nc>%1.1f)', ...
                        subjnrs(sj), whichModelToPlot, ve_thresh, nc_thresh))
                else
                    sgtitle(sprintf('Subject S%d: Seq vs Sim amplitudes (voxels with pRF ve>%1.1f & nc>%1.1f)', ...
                        subjnrs(sj),ve_thresh, nc_thresh))
                end
            else
                if plotModelAmpl
                    sgtitle(sprintf('Subject S%d: Seq vs Sim predicted amplitudes (all voxels)',subjnrs(sj)))
                else
                    sgtitle(sprintf('Subject S%d: Seq vs Sim amplitudes (all voxels)',subjnrs(sj)))
                end
            end
            
            % Save figure
            if saveFigs
                if plotModelAmpl
                    fName = sprintf('Summary%sSubject%d_SEQvsSIM_%sModelAmpl_XVal1_%s',...
                        replace(string(roisToPlot(idx)),'/',''), subjnrs(sj), whichModelToPlot, roiType);
                else
                    fName = sprintf('Summary%sSubject%d_SEQvsSIMAmpl_XVal1_%s',...
                        replace(string(unique(tmpT.ROI)),'/',''), subjnrs(sj),roiType);
                end
                
                if plotFitFlag
                    fName = sprintf('%s_LMMfit_%s',fName, LMMlbl);
                    if plotDataFlag
                        fName = [fName '_wData'];
                    else
                        fName = [fName '_noData'];
                    end
                else
                    fName = sprintf('%s_dataOnly',fName);
                end
                thisSaveFigDir = fullfile(saveFigDir);
                if ~exist(thisSaveFigDir,'dir'); mkdir(thisSaveFigDir); end
                saveas(gcf, fullfile(thisSaveFigDir, [fName '.png']))
%                 print(gcf,fullfile(thisSaveFigDir,fName),'-depsc')
            end
        end % ROIS
    end % SUBJNRS
end % PLOTSUBJTOGETHER FLAG