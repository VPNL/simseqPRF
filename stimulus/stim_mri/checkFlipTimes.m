function [] = checkFlipTimes(fileNameToLoad)
% Script to ad hoc check stimulus frame timing.

% Example:
% run = 999;
% date = '20201125T122751';
% fName = sprintf('simseqPRF_test_run%d_%s.mat',run, date);
% checkFlipTimes(fName)
set(0,'DefaultFigureWindowStyle','docked') % use 'normal' to undock

% Load output file 
load(fullfile(fileNameToLoad));

desiredFlipTime = diff(scan.allFlips);
measuredFlipTime  = diff(response.VBLTimestamp); % same as response.stimOnsetTime

% Set y range 
yrange = 0.02 .* [-1 1];

% Plot!
figure(101); clf
plot(desiredFlipTime, 'g', 'LineWidth',3);
hold on; plot(measuredFlipTime, 'r-','LineWidth',2);

% Make plot pretty
xlabel('flip #')
ylabel('Measured flip time (s)')
set(gca, 'TickDir', 'out', 'FontSize', 20)
ylim(median(measuredFlipTime) + yrange)
xlim([0 length(desiredFlipTime)])
legend('Desired flip time', 'Measured flip time'); legend box off

% Print frames between stimuli
frames = round(diff(measuredFlipTime) / (1/60));

% how many interstimulus frames differed from the median?
fprintf('Nr of off-time flips: \t %3.0f\n', sum(frames ~= median(frames))/2)