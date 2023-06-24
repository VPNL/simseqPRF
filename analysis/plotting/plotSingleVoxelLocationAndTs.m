function fH = plotSingleVoxelLocationAndTs(roiVoxelData, prfData, onsetTimeSecs, selectedVoxels, TR, conditions, pths, roiName, subFolder, saveFigs, savePth)
%% Plot location and timseries of single voxel 
fH = figure(99); clf; set(gcf, 'position', [73,555,1904,555], 'color', 'w')

roisplit = strsplit(roiName, '_');
roi = roisplit{2}; 
if strcmp(roi, 'lh')
    roi = roisplit{3};
    clear roisplit;
end

% stim coords
loc = getStimLocDeg(pths.projectDir, pths.subjnr);
f = [1 2 3 4];
for vv = f
    v{vv} = [loc.xBounds(vv+1,1) loc.yBounds(vv+1,1); ...
        loc.xBounds(vv+1,2) loc.yBounds(vv+1,1);
        loc.xBounds(vv+1,2) loc.yBounds(vv+1,2);
        loc.xBounds(vv+1,1) loc.yBounds(vv+1,2)];
end

% Get time axis
t = 0:(size(roiVoxelData,1)-1);

% Define onset of conditions
tx = [onsetTimeSecs size(roiVoxelData,1)*TR];
condColors = {[1 1 1],[0.5450 0 0],[0 0 0.5450]};

meanRoiVoxelData = squeeze(mean(roiVoxelData,2, 'omitnan'));

% Plot selected voxels
for vv = selectedVoxels
    subplot(1,2,1); cla 
    gx = gca;
    gxPos = get(gx,'Position'); set(gx,'Units','normalized')
    set(gx,'Position',[gxPos(1) gxPos(2)*1.15 gxPos(3) gxPos(4)])
    
    axis equal;  xlim([-1,12]); ylim([-1,12]);
    set(gx, 'FontSize', 15, 'TickDir', 'out', ...
        'XTick',[0 5 10], 'XTickLabel',{'0','5','10'}, ...
        'YTick',[0 5 10], 'YTickLabel',{'0','5','10'});
    xlabel({'x position (deg)',''});
    ylabel('y position (deg)');
    title(sprintf('%s: voxel %d, prf ve = %1.2f',roi, vv, prfData.varexpl(vv)))
    
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
    mx = max(abs([min(meanRoiVoxelData(:,vv)),max(meanRoiVoxelData(:,vv))]));
    if mx > 1
        yl = [-1,1] .* (1.1*mx);
    else
        yl = [-1, 1];
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
        ind = conditions(ti);
        C = condColors{ind+1};
        patch(X,Y,C,'FaceAlpha', 0.5, 'EdgeAlpha',0); hold on;
    end
    
    % Add individual voxel time series and mean
    plot(t,zeros(1,length(t)),'k'); hold on;
    plot(t,meanRoiVoxelData(:,vv), 'color', 'k', 'LineWidth',3)
    
    % Make plot pretty
    set(gca, 'FontSize', 15, 'TickDir', 'out')
    
    ylabel('% Signal change');
    xlabel('TRs (s)')
    title(sprintf('%s: average time series', roi))
    
    if saveFigs
        if ~exist(fullfile(savePth, sprintf('voxelTs_%s',roiName),subFolder),'dir')
            mkdir(fullfile(savePth, sprintf('voxelTs_%s',roiName),subFolder));
        end
        saveas(fH,fullfile(savePth, sprintf('voxelTs_%s',roiName),subFolder, sprintf('voxel%d_ts.png',vv)))
    end
end
