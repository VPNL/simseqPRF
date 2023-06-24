  %% runme_simSeqPRF
% Add toolboxes: tbUse({'Psychtoolbox-3', 'simseqPRF'})
% at the scanner:
% set trigger box to these options (per JG):
% change modes / manual configure / HHSC 1X5D / USB / HID NAR NO5

subject       = 'test';
uniqueRuns    = 2;
totalRuns     = 2;
versionNr     = 5;      % each subject has its own version nr
offset        = [0 0];  % center offset of image [X,Y] in pixels
atScanner     = 0;      % 0 = don't start sequence with "t", 1 = do start sequence with "t"
coil          = 16;     % 16 or 32-chan coil, 16 will not flip up/down and flip L/R
useEyetracker = 0;
loadTrials    = [];     % load the specific stimuli and order from trialSets

if ~exist(fullfile(simseqRootPath, 'data'),'dir'); mkdir(fullfile(simseqRootPath, 'data')); end

for r = 1:totalRuns

    runNum = mod(r-1,uniqueRuns)+1;
   
    singleRun_simSeqPRF_varySize_x_Duration_mainExp(subject, ...
        runNum,...
        versionNr, ...
        loadTrials,...
        offset,...
        atScanner,...
        coil, ...
        useEyetracker)

    % Check flip times
    dateToday = datestr(now,'yyyymmdd');
    d = dir(fullfile(simseqRootPath,'data',['S' subject],sprintf('simseqPRF_%s_run%d_v%d_%s*.mat',subject, r, versionNr, dateToday)));
    if ~isempty(d)
        checkFlipTimes(fullfile(d(1).folder, d(1).name))
    end
    clear d
end