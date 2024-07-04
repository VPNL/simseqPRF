function fH = makeSupplFigure5_eyefixation(projectDir,subjnrs,saveFigs,saveFigDir)
% Function to reproduce supplementary figure 5: 
% Plotting the distribution of eye gaze position
%
% From the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University
% 
% INPUTS (required):
% - projectDir      : base folder of project
% - subjnrs         : Subjects to plot
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%% Check inputs
if ~exist('subjnrs','var') || isempty(subjnrs)
    subjnrs = [1,2,3,4,9];
end
if ~ismember(subjnrs,[1,2,3,7,12])
    error('[%]: Check subject nrs, not all subjects have eye tracking data',mfilename)
end
    
%% Visualize fixation points from subject's eyetracking data
fprintf('Loading subject data from file\n')
pths = getSubjectPaths(projectDir,1);
loadDataDir = fullfile(simseqRootPath,'data','eyetracking');
load(fullfile(loadDataDir,[sprintf('S%i_',subjnrs) 'eyetrackingDataClean.mat']))

if saveFigs
    saveFigDir = fullfile(simseqRootPath,'results', 'group', ...
        'eyetrackingPerformance_SimSeq');
    if ~exist(saveFigDir, 'dir'); mkdir(saveFigDir); end
end

% check percentage of data
prctData = NaN(size(allSubjectData,1),size(allSubjectData,2));
for sj = 1:size(allSubjectData,1)
    for r = 1:size(allSubjectData,2)
        if ~isempty(allSubjectData(sj,r).xyDataDeg)
            prctData(sj,r) = 1-sum(isnan(allSubjectData(sj,r).xyDataDeg(:,1)))/length(allSubjectData(sj,r).xyDataDeg(:,1));
        end
    end
end
% Keep runs with more than 20% data.
toKeep = (prctData>0.2);

%% Compute group mean, median and covariance
subjXYPosMean   = NaN(size(allSubjectData,1),size(allSubjectData,2),2);
subjXYPosMedian = subjXYPosMean;
subjXYPosCov    = NaN(size(allSubjectData,1),size(allSubjectData,2),2,2);
allSubjPos      = [];

for sj = 1:size(allSubjectData,1)
     for r = 1:size(allSubjectData,2)
       if toKeep(sj,r)==1
            subjXYPosMean(sj,r,:) = mean(allSubjectData(sj,r).xyDataDeg,1,'omitnan');
            subjXYPosMedian(sj,r,:) = median(allSubjectData(sj,r).xyDataDeg,1,'omitnan');
            subjXYPosCov(sj,r,:,:) = cov(allSubjectData(sj,r).xyDataDeg(~isnan(allSubjectData(sj,r).xyDataDeg(:,1)),:));
            allSubjPos = cat(1,allSubjPos, allSubjectData(sj,r).xyDataDeg);
       end
     end
end
    
mnCovRun = reshape(squeeze(mean(subjXYPosCov,2, 'omitnan')),5,2,2);
mnMdRun = squeeze(mean(subjXYPosMedian,2,'omitnan'));

%% Plot kernel density of fixations (individual subjects + group)
fH = figure; set(gcf,'Color', 'w','Position',[1 1 1600 1200]); 

for s = 1:size(allSubjectData,1)
    
    subjData = [];
    for r = 1:8
        subjData = cat(1,subjData, allSubjectData(s,r).xyDataDeg);
    end
    % distribution of samples
    [~,density,X,Y] = kde2d(subjData,2^10, [-5,-5],[5 5]);
    
    subplot(2,3,s); cla
    
    % Plot Stim
    p1 = patch([0.59, 0.59, 4.59, 4.59],[0.59, 4.59, 4.59, 0.59],'k'); hold all;
    p1.FaceAlpha = 0.1;
    p2 = patch([2.59, 2.59, 4.59, 4.59],[2.59, 4.59, 4.59, 2.59],'k');
    p2.FaceAlpha = 0.2;
    p3 = patch(-1.*[0.59, 0.59, 4.59, 4.59],-1.*[0.59, 4.59, 4.59, 0.59],'k');
    p3.FaceAlpha = 0.1;
    p4 = patch(-1.*[2.59, 2.59, 4.59, 4.59],-1.*[2.59, 4.59, 4.59, 2.59],'k');
    p4.FaceAlpha = 0.2;

    % Plot contours
    contourf(X,Y,density./max(density(:)),[0.01,0.1, 0.5,1],'k'); hold on;
    c = colorbar; c.Ticks = [0.01,0.1,0.5,1];

    % Can also be surface
    %     surf(X,Y,density./max(density(:)),'EdgeColor','none'); hold all;
    %     view(2)    
    % Or plot 95% confidence ellipse
    %     error_ellipse(squeeze(mnCovRun(s,:,:)),squeeze(mnMdRun(s,:)),'conf',0.95,'color','k'); hold on;
    
    % Plot median as red cross
    plot(mnMdRun(s,1),mnMdRun(s,2),'+r','markersize',2); hold on;
    title(sprintf('S%02d',subjnrs(s)));
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    
    if ismember(s,[4,5]),xlabel('Horizontal (deg)'); end
    if ismember(s, [1,4]), ylabel('Vertical (deg)'); end
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
end   
   
% Now add Group
% Squeeze nans out of matrix to save some memory
x0All = allSubjPos(~isnan(allSubjPos(:,1)),1);
y0All = allSubjPos(~isnan(allSubjPos(:,1)),2);

% Compute kernel density
[~,density,X,Y] = kde2d([x0All,y0All],2^10, [-5 -5],[5 5]);

% Plot Stim
subplot(2,3,6)
p1 = patch([0.59, 0.59, 4.59, 4.59],[0.59, 4.59, 4.59, 0.59],'k'); hold all;
p1.FaceAlpha = 0.1;
p2 = patch([2.59, 2.59, 4.59, 4.59],[2.59, 4.59, 4.59, 2.59],'k');
p2.FaceAlpha = 0.2;
p3 = patch(-1.*[0.59, 0.59, 4.59, 4.59],-1.*[0.59, 4.59, 4.59, 0.59],'k');
p3.FaceAlpha = 0.1;
p4 = patch(-1.*[2.59, 2.59, 4.59, 4.59],-1.*[2.59, 4.59, 4.59, 2.59],'k');
p4.FaceAlpha = 0.2;

% Plot all subject data (N=5)
contourf(X,Y,density./max(density(:)),[0.01,0.1,0.5,1]);
c=colorbar; c.Ticks = [0.01,0.5,1];

% Can also be a surface
% surf(X,Y,density./max(density(:)),'EdgeColor','none'); view(2);

grid on; axis square;
xlim(5*[-1,1]); ylim(5*[-1,1]);
xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
title('N=5')

% Plot median
plot(mean(mnMdRun(:,1),1,'omitnan'),mean(mnMdRun(:,2),1,'omitnan'),'+r','markersize',5); hold on;

if saveFigs
    fName = sprintf('supplFig5_eyetracking_XYDensity_wContours');
    if ~exist(fullfile(saveFigDir,'supplFig5'),'dir')
        mkdir(fullfile(saveFigDir,'supplFig5'));
    end
    saveas(gcf,fullfile(saveFigDir,'supplFig5',[fName '.png']))
    %     print(gcf,fullfile(saveFigDir,fName),'-depsc','-painters','-r300','-loose');
end

return

