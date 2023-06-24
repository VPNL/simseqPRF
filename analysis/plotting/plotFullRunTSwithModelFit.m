function fH = plotFullRunTSwithModelFit(fH, dataToPlot, modelPredToPlot, ...
    stimToPlot, ttl)

% Check inputs
if ~exist('fH', 'var') || isempty(fH)
    fH = figure(101); clf; set(gcf, 'position', [55,440,1994,522]);
end

if ~exist('stimToPlot','var') || isempty(stimToPlot)
    stimToPlot = [];
end

% Set up figure and axes
figure(fH); hold on; clf;

yl = 1.1.*[min(dataToPlot),max(dataToPlot)];
if isnan(yl); yl = [1 1]; end

vshift = yl(2) + 0.1;
yl(2) = vshift + 1.1;

numTimePoints = size(dataToPlot,1);
t = 0:1:(numTimePoints-1);
tInMs = 0:0.001:numTimePoints; 
tInMs = tInMs(1:end-1);
makeprettyfigures;

% Plot it!
plot(t,zeros(1,length(t)),'k', 'LineWidth',0.5); hold on;
if ~isempty(stimToPlot)
    plot(tInMs,stimToPlot+vshift, 'color', [0.5 0.5 0.5 0.5], 'LineWidth',1);
end
plot(t,dataToPlot, 'ko-');
plot(t,modelPredToPlot,'color','r','lineWidth',2.5);
ylim(yl); 
xlim([0 numTimePoints]); 
box off;

title(ttl, 'Interpreter','none')

l = findobj('Type','line');
legend(l([2,1]),{'Average Data','Model prediction'}, 'Location', 'SouthEast'); legend boxoff;
xlabel('Time (TRs)'); ylabel('BOLD signal (% change)');

return
