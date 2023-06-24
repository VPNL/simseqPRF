function fH = plotVoxelpRFLocationDataAndModelTs(roiVoxelData, prfData, modelPredData, varargin)
% fH = plotVoxelpRFLocationAndTs(roiVoxelData, prfData, modelPredData, ...
%           [TR], [conditionData], [roiName], [saveFigs], [savePth])
%
% Function to plot a single vertex or vertices pRF location next to it's time series
%
% INPUTS
% roiVoxelData          : matrix with time series [time x voxels]
% prfData               : struct with prf properties as fields
%                           - xData
%                           - yData
%                           - size (or sigma). Note: function does not correct for cssFit size changes
%                           - variance explained
% modelPredData         : struct with model prediction properties
%                           - ts (predicted time course)
%                           - beta (beta values for temporal channels)
%                           - R2
%                           - temporalModel (name of used temporal model)
%                           - spatialModel (name of used temporal model)
% [TR]                  : integer with TR in seconds (default = 1s)
% [conditionData]       : struct with condition number, onset times, colors
% [roiName]             : string defining ROI name
% [saveFigs]            : boolean to save figures or not
% [savePth]             : where to save figures
%
% OUTPUTS
% fH                    : Figure handle

%% Parse inputs
p = inputParser;
p.addRequired('roiVoxelData', @isnumeric);
p.addRequired('prfData',@isstruct);
p.addRequired('modelPredData',@isstruct);
p.addParameter('TR',1, @isnumeric);
p.addParameter('conditionData',[], @isstruct);
p.addParameter('roiName', 'ROI', @ischar);
p.addParameter('plotAvg', false, @islogical);
p.addParameter('saveFigs', false, @islogical);
p.addParameter('savePth', pwd, @ischar);

p.parse(roiVoxelData,prfData,modelPredData,varargin{:});

% Rename variables
roiVoxelData   =  p.Results.roiVoxelData;
prfData        =  p.Results.prfData;
modelPredData  =  p.Results.modelPredData;
TR             =  p.Results.TR;
conditionData  =  p.Results.conditionData;
roiName        =  p.Results.roiName;
plotAvg        =  p.Results.plotAvg;
saveFigs       =  p.Results.saveFigs;
savePth        =  p.Results.savePth;

if ~isempty(conditionData)
    % Get stim coords
    f = [1 2 3 4]; % define faces
    for vv = f % define vertices
        v{vv} = [conditionData.loc.xBounds(vv+1,1) conditionData.loc.yBounds(vv+1,1); ...
            conditionData.loc.xBounds(vv+1,2) conditionData.loc.yBounds(vv+1,1);
            conditionData.loc.xBounds(vv+1,2) conditionData.loc.yBounds(vv+1,2);
            conditionData.loc.xBounds(vv+1,1) conditionData.loc.yBounds(vv+1,2)];
    end
    % Define onset of conditions
    tx = [conditionData.onsetTimeSecs size(roiVoxelData,1)*TR];
end

% Define colors for beta value bars
betaColors = [0.3 0.3 0.3; 0.7 0.7 0.7];

if plotAvg
    % Get mean voxel data
    dataToPlot       = squeeze(mean(roiVoxelData,2, 'omitnan'));
    predictionToPlot = squeeze(mean(modelPredData.ts,2, 'omitnan'));
    betasToPlot      = squeeze(mean(modelPredData.beta,2, 'omitnan'));
else
    dataToPlot       = roiVoxelData;
    predictionToPlot = modelPredData.ts;
    betasToPlot      = modelPredData.beta;
end

% Get time axis (TRs)
t = 0:TR:(size(dataToPlot,1)-1);

% Plot location and timseries of single voxel
fH = figure(99); clf; set(gcf, 'position', [73,555,1904,555], 'color', 'w')

% Plot selected voxels
if plotAvg
    
    subplot(1,4,1); hold all;
    gx = gca;
    gxPos = get(gx,'Position'); set(gx,'Units','normalized')
    set(gx,'Position',[gxPos(1) gxPos(2)*1.15 gxPos(3) gxPos(4)])
    
    axis equal;  xlim([-1,12]); ylim([-1,12]);
    set(gx, 'FontSize', 15, 'TickDir', 'out', ...
        'XTick',[0 5 10], 'XTickLabel',{'0','5','10'}, ...
        'YTick',[0 5 10], 'YTickLabel',{'0','5','10'});
    xlabel({'x position (deg)',''});
    ylabel('y position (deg)');
    title(sprintf('All pRFs'))
    
    % Plot stimulus square
    for jj = 1:length(v)
        p = patch('Faces',f,'Vertices',v{jj},'FaceColor','red'); hold on;
    end
    
    gx.Units = 'Points';
    gxSizePix = gx.Position(3:4);
    gx.Units = 'normalized';
    pixPerUnit = ceil(max(gxSizePix ./ [range(xlim(gx)),range(ylim(gx))]));
    
    for vv = 1:size(roiVoxelData,2)
        % Plot pRF
        plot([-1,12], [0,0],'k', [0,0],[-1,12], 'k');
        scatter(prfData.x0(vv),prfData.y0(vv),2*(prfData.size(vv)*pixPerUnit).^2, ...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[1 1 .5],...
            'MarkerFaceAlpha',0.2, ...
            'LineWidth',2)
    end
    
    %% reset ylim
    yl_range = [min(dataToPlot),max(dataToPlot)];
    yl = [min(yl_range)-0.1, max(yl_range)+0.1];
    if sum(yl)==0 || yl(1) > yl(2) || any(isnan(yl))
        yl = [-1 1];
    end
    
    % Plot time series + condition markers
    subplot(1,4,[2 3]); cla;
    ylim(yl); xlim([min(t),max(t)]);
    AX = axis;
    for ti = 1:length(tx)-1
        t1 = tx(ti);
        t2 = tx(ti+1);
        X = [t1 t1 t2 t2];
        Y = [AX(3) AX(4) AX(4) AX(3)];
        ind = conditionData.conditions(ti);
        C = conditionData.condColors{ind+1};
        patch(X,Y,C,'FaceAlpha', 0.5, 'EdgeAlpha',0); hold on;
    end
    
    % Add individual voxel time series and mean
    plot(t,zeros(1,length(t)),'k', 'lineWidth',0.5); hold on;
    plot(t,dataToPlot, 'ko-', 'LineWidth',1.5)
    plot(t,predictionToPlot, 'r', 'lineWidth',3)
    
    % Make plot pretty
    set(gca, 'FontSize', 15, 'TickDir', 'out')
    
    ylabel('% Signal change');
    xlabel('TRs (s)')
    title(sprintf('%s: %s %s R^2: %2.1f', roiName, ...
            modelPredData.spatialModel, modelPredData.temporalModel, ...
            100*mean(modelPredData.R2, 'omitnan')))
    
    subplot(1,4,4); cla; 
    ax = gca;
    set(ax,'Position',[ax.Position(1), ax.Position(2), ax.Position(3)*0.8, ax.Position(4)*0.8])
    hold on;
    if strcmp(modelPredData.temporalModel,'2ch-exp-sig')
        xlabels = {'S-exp','T-sig'};
        xl = [0.5 2.5];
        bb = [1,2];
    elseif strcmp(modelPredData.temporalModel,'2ch-css-sig')
        xlabels = {'S-css','T-sig'};
        xl = [0.5 2.5];
        bb = [1,2];
    else
        xlabels = {'1ch linear'};
        xl = [0.5 1.5];
        bb = 1;
    end

    % if we want each bar to have a different color, loop
    for b = bb
        bar(b, betasToPlot(b), 'FaceColor',  betaColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);
    end

    % show standard deviation on top
    semsB = std(modelPredData.beta,2, 'omitnan')./sqrt(size(dataToPlot,2));
    h = ploterr(bb, betasToPlot, [], semsB, 'k.', 'abshhxy', 0);
    set(h(1), 'marker', 'none'); % remove marker

    yl = 1.1*[min(betasToPlot(:)-semsB(:)),max(betasToPlot(:)+semsB(:))];
    if yl(1)>0; yl(1) = 0; end
    if yl(2)<0; yl(2) = 0; end
    set(gca, 'xtick', bb, 'xticklabel', xlabels, 'xlim', xl, 'ylim',yl);
    ylabel('Beta weight'); box off;
    title({sprintf('%s',roiName),sprintf('%s',modelPredData.spatialModel)})
    
    
    if saveFigs
        if ~exist(fullfile(savePth, sprintf('voxelTs_%s',roiName)),'dir')
            mkdir(fullfile(savePth, sprintf('voxelTs_%s',roiName)));
        end
        plotName = 'voxelAvg_ts.png';
        saveas(fH,fullfile(savePth, sprintf('voxelTs_%s',roiName), plotName))
    end
    
    
else
    if ~exist('prfData.exponent','var')
        prfData.exponent = NaN(size(prfData.x0));
    end
    for vv = 1:size(dataToPlot,2)
        
        subplot(1,4,1); cla;
        gx = gca;
        gxPos = get(gx,'Position'); set(gx,'Units','normalized')
        set(gx,'Position',[gxPos(1) gxPos(2)*1.15 gxPos(3) gxPos(4)])
        
        axis equal;  xlim([-1,12]); ylim([-1,12]);
        set(gx, 'FontSize', 15, 'TickDir', 'out', ...
            'XTick',[0 5 10], 'XTickLabel',{'0','5','10'}, ...
            'YTick',[0 5 10], 'YTickLabel',{'0','5','10'});
        xlabel({'x position (deg)',''});
        ylabel('y position (deg)');
        
        title(sprintf('pRF %d, ve=%2.1f%%, sd=%1.1f deg, exp=%1.2f', vv, prfData.varexpl(vv)*100,prfData.size(vv),prfData.exponent(vv)))
        
        % Plot stimulus square
        for jj = 1:length(v)
            p = patch('Faces',f,'Vertices',v{jj},'FaceColor','red'); hold on;
        end
        
        gx.Units = 'Points';
        gxSizePix = gx.Position(3:4);
        gx.Units = 'normalized';
        pixPerUnit = ceil(max(gxSizePix ./ [range(xlim(gx)),range(ylim(gx))]));
        
        % Plot pRF
        plot([-1,12], [0,0],'k','LineWidth',0.5); hold on
        plot([0,0],[-1,12], 'k', 'LineWidth',0.5);
        if prfData.size(vv) > 10.5
            xlim([-1,20]); ylim([-1,20]);
        end
        scatter(prfData.x0(vv),prfData.y0(vv),2*(prfData.size(vv)*pixPerUnit).^2, ...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[1 1 .5],...
            'MarkerFaceAlpha',0.2, ...
            'LineWidth',2)
        
        %% reset ylim
        yl_range = ([min(dataToPlot(:,vv)),max(dataToPlot(:,vv))]);
        yl = [min(yl_range)-0.1, max(yl_range)+0.1];
        if sum(yl)==0
            yl = [-1 1];
        end
        % Plot time series + condition markers
        subplot(1,4,[2 3]); cla;
        ax = gca;
        set(ax,'Position',[ax.Position(1), ax.Position(2)+0.05, ax.Position(3), ax.Position(4)-0.05])
        ylim(yl); xlim([min(t),max(t)]);
        AX = axis;
        for ti = 1:length(tx)-1
            t1 = tx(ti);
            t2 = tx(ti+1);
            X = [t1 t1 t2 t2];
            Y = [AX(3) AX(4) AX(4) AX(3)];
            ind = conditionData.conditions(ti);
            C = conditionData.condColors{ind+1};
            patch(X,Y,C,'FaceAlpha', 0.5, 'EdgeAlpha',0); hold on;
        end
        
        % Add individual voxel time series and mean
        plot(t,zeros(1,length(t)),'k', 'LineWidth',0.5); hold on;
        plot(t,dataToPlot(:,vv), 'ko-', 'LineWidth',1.5)
        plot(t,predictionToPlot(:,vv), 'r', 'LineWidth',3)
        % Make plot pretty
        set(gca, 'FontSize', 15, 'TickDir', 'out')
        
        ylabel('% Signal change');
        xlabel('TRs (s)')
        title(sprintf('%s %s %s R^2=%2.1f%%', roiName, ...
            modelPredData.spatialModel, modelPredData.temporalModel, ...
            modelPredData.R2(vv)*100))
        
        subplot(1,4,4); cla; 
        ax = gca;
        set(ax,'Position',[ax.Position(1)+0.01, ax.Position(2), ax.Position(3)*0.8, ax.Position(4)*0.8])
        hold on;
        if strcmp(modelPredData.temporalModel,'2ch-exp-sig')
            
            beta = betasToPlot(vv,:);
            bb = 1:length(beta);
            xl = [0.5 length(beta)+0.5];
            if bb == 1
                xlabels = {'S+T'};
            else
                xlabels = {'S-exp','T-sig'};
            end
        elseif strcmp(modelPredData.temporalModel,'2ch-css-sig')
            
            beta = betasToPlot(vv,:);
            bb = 1:length(beta);
            xl = [0.5 length(beta)+0.5];
            if bb == 1
                xlabels = {'S+T'};
            else
                xlabels = {'S-css','T-sig'};
            end
        else
            xlabels = {'1ch linear'};
            xl = [0.5 1.5];
            bb = 1;
            beta = betasToPlot(vv);
        end
        
        % if we want each bar to have a different color, loop
        for b = bb
            bar(b, beta(b), 'FaceColor',  betaColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);
        end
        
        yl = [-1, 1]; %1.1*[min(beta(:)),max(beta(:))];
         if min(beta(:))<yl(1), yl(1) = 1.1*min(beta(:)); end
         if max(beta(:))>yl(2), yl(2) = 1.1*max(beta(:)); end
%          if isnan(yl),yl = [0,1]; end
        set(gca, 'xtick', bb, 'xticklabel', xlabels, 'xlim', xl, 'ylim',yl);
        ylabel('Beta weight'); box off;
        
        if saveFigs
            if ~exist(fullfile(savePth, sprintf('voxelTs_%s',roiName)),'dir')
                mkdir(fullfile(savePth, sprintf('voxelTs_%s',roiName)));
            end
            plotName = sprintf('voxel%d_ts.png',vv);
            saveas(fH,fullfile(savePth, sprintf('voxelTs_%s',roiName), plotName))
        end
    end
end
return
