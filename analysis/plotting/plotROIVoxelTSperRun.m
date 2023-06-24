function fH = plotROIVoxelTSperRun(roiName, roiVoxelData, onsetTimeSecs, TR, conditions, saveFigs, savePth)
%% Plot all ROI voxel time series, separate per run

fH = figure; clf; set(gcf, 'Position', [704,52,1237,925], 'color', 'w'); hold all;

% Get time axis
t = 0:(size(roiVoxelData,1)-1);

% Define onset of conditions
tx = [onsetTimeSecs size(roiVoxelData,1)*TR];

yl = [-10,10];
condColors = {[1 1 1],[0.5450 0 0],[0 0 0.5450]};
cmap = lines(size(roiVoxelData,2));

for r = 1:size(roiVoxelData,2)
    
    subplot(size(roiVoxelData,2),1,r); hold all
    
    % Make plot pretty
    set(gca, 'FontSize', 15, 'TickDir', 'out')
    xlim([min(t),max(t)]); ylim(yl);
    AX = axis;
    title(sprintf('Scan %d',r))
    if r == 3; ylabel('% Signal change'); end
    
    % make patches for the different trial conditions
    % (adapted from tc_plotWholeTc.m)
    for ii = 1:length(tx)-1
        t1 = tx(ii);
        t2 = tx(ii+1);
        X = [t1 t1 t2 t2];
        Y = [AX(3) AX(4) AX(4) AX(3)];
        ind = conditions(ii);
        C = condColors{ind+1};
        patch(X,Y,C,'FaceAlpha', 0.5, 'EdgeAlpha',0);
    end
    
    % Add individual voxel time series and mean
    plot(t,squeeze(roiVoxelData(:,r,:)), 'color', cmap(r,:), 'LineWidth',0.75)
    plot(t,mean(roiVoxelData(:,r,:),3), 'color', 'k', 'LineWidth',3)
    
end
% Only plot xlabel at the bottom to keep figure labels sparse
xlabel('TRs (s)')


if saveFigs
    print(fullfile(savePth, sprintf('voxelTs_%s',roiName),...
        sprintf('allVoxelsTsPerScan')),'-dpng')
end