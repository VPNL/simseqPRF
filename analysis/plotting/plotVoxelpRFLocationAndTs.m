function fH = plotVoxelpRFLocationAndTs(roiVoxelData, prfData, varargin)
% fH = plotVoxelpRFLocationAndTs(roiVoxelData, prfData, [TR], [conditionData], [roiName], [saveFigs], [savePth])
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
p.addParameter('TR',1, @isnumeric);
p.addParameter('conditionData',[], @isstruct);
p.addParameter('roiName', 'ROI', @ischar);
p.addParameter('plotAvg', true, @islogical);
p.addParameter('saveFigs', false, @islogical);
p.addParameter('savePth', pwd, @ischar);

p.parse(roiVoxelData,prfData,varargin{:});

% Rename variables
roiVoxelData   =  p.Results.roiVoxelData;
prfData        =  p.Results.prfData;
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


if plotAvg
    % Get mean voxel data
    dataToPlot = squeeze(mean(roiVoxelData,2, 'omitnan'));
else
    dataToPlot = roiVoxelData;
end

% Get time axis (TRs)
t = 0:TR:(size(dataToPlot,1)-1);

% Plot location and timseries of single voxel
fH = figure(99); clf; set(gcf, 'position', [73,555,1904,555], 'color', 'w')

% Plot selected voxels
if plotAvg
    
    subplot(1,2,1); hold all;
    gx = gca;
    gxPos = get(gx,'Position'); set(gx,'Units','normalized')
    set(gx,'Position',[gxPos(1) gxPos(2)*1.15 gxPos(3) gxPos(4)])
    
    axis equal;  xlim([-1,12]); ylim([-1,12]);
    set(gx, 'FontSize', 15, 'TickDir', 'out', ...
        'XTick',[0 5 10], 'XTickLabel',{'0','5','10'}, ...
        'YTick',[0 5 10], 'YTickLabel',{'0','5','10'});
    xlabel({'x position (deg)',''});
    ylabel('y position (deg)');
    title(sprintf('%s: all pRFs', roiName))
    
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
    subplot(1,2,2); cla;
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
    plot(t,zeros(1,length(t)),'k'); hold on;
    plot(t,dataToPlot, 'color', 'k', 'LineWidth',3)
    
    % Make plot pretty
    set(gca, 'FontSize', 15, 'TickDir', 'out')
    
    ylabel('% Signal change');
    xlabel('TRs (s)')
    title(sprintf('%s: average time series', roiName))
    
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
        
        subplot(1,2,1); cla;
        gx = gca;
        gxPos = get(gx,'Position'); set(gx,'Units','normalized')
        set(gx,'Position',[gxPos(1) gxPos(2)*1.15 gxPos(3) gxPos(4)])
        
        axis equal;  xlim([-1,12]); ylim([-1,12]);
        set(gx, 'FontSize', 15, 'TickDir', 'out', ...
            'XTick',[0 5 10], 'XTickLabel',{'0','5','10'}, ...
            'YTick',[0 5 10], 'YTickLabel',{'0','5','10'});
        xlabel({'x position (deg)',''});
        ylabel('y position (deg)');
        
        title(sprintf('%s: voxel %d, prf ve = %1.2f, size = %1.2f, exp = %1.2f',roiName, vv, prfData.varexpl(vv),prfData.size(vv),prfData.exponent(vv)))
        
        % Plot stimulus square
        for jj = 1:length(v)
            p = patch('Faces',f,'Vertices',v{jj},'FaceColor','red'); hold on;
        end
        
        gx.Units = 'Points';
        gxSizePix = gx.Position(3:4);
        gx.Units = 'normalized';
        pixPerUnit = ceil(max(gxSizePix ./ [range(xlim(gx)),range(ylim(gx))]));
        
        % Plot pRF
        plot([-1,12], [0,0],'k', [0,0],[-1,12], 'k');
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
        subplot(1,2,2); cla;
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
        plot(t,zeros(1,length(t)),'k'); hold on;
        plot(t,dataToPlot(:,vv), 'color', 'k', 'LineWidth',3)
        
        % Make plot pretty
        set(gca, 'FontSize', 15, 'TickDir', 'out')
        
        ylabel('% Signal change');
        xlabel('TRs (s)')
        title(sprintf('%s: average time series', roiName))
        
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
