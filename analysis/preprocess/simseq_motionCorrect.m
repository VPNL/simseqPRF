function simseq_motionCorrect(sessionDir)
% simseq_motionCorrect(sessionDir)
% 
% Function to apply within and between scan motion correction for simseq experiment
%
% Within scan motion correction is relative to the 8th frame of each scan
% (see params.motionCompRefFrame = 8; % reference TR for motion correction)
% Between scan motion correction relative to the first scan in the session
% (see baseScan = 1);
%
%
% Inputs:
% sessionDir    : (str) path to subject's VistaSoft session
%
% Written by EK 2021, adapted from D Finzi toonotopy code

cd(sessionDir);
fns = dir(sessionDir);
fns = {fns.name};
lid = fopen('preprocess_log.txt', 'a+');
fprintf(lid, 'Starting motion analysis for session %s. \n\n', sessionDir);

%% Load initialized session
load mrInit_params.mat
load mrSESSION.mat

% Open hidden inplane
% hi = initHiddenInplane('Original', 1);
hi = initHiddenInplane('Timed', 1);

% Update params
params.motionCompRefFrame = 8; % reference TR for motion correction
params.motionCompSmoothFrames = 3; % smoothing window for motion correction

% get number of runs
[tmp runCount]=size(dataTYPES(1).scanParams);

% do within-scan motion compensation (unless its already been done)
fprintf(lid, 'Starting within-scan motion compensation... \n');
fprintf('Starting within-scan motion compensation... \n');
setpref('VISTA', 'verbose', false); % suppress wait bar
if ~exist(fullfile(sessionDir, 'Images', 'Within_Scan_Motion_Est.fig'), 'file')
    hi = motionCompSelScan(hi, 'MotionComp', 1:runCount, ...
        params.motionCompRefFrame, params.motionCompSmoothFrames);
    saveSession; % close all;
else 
    warning('Looks within-scan motion was already completed. Will skip this step')
end

fprintf(lid, 'Within-scan motion compensation complete. \n\n');
fprintf('Within-scan motion compensation complete. \n\n');

% do between-scan motion compensation (unless its already been done)
fprintf(lid, 'Starting between-scan motion compensation... \n');
fprintf('Starting between-scan motion compensation... \n');
if ~exist(fullfile(sessionDir, 'Between_Scan_Motion.txt'), 'file')
    hi = initHiddenInplane('MotionComp', 1);
    baseScan = 1;
    targetScans = 1:runCount;
    [hi, M] = betweenScanMotComp(hi, 'MotionComp_RefScan1', baseScan, targetScans);
    fname = fullfile('Inplane', 'MotionComp_RefScan1', 'ScanMotionCompParams');
    save(fname, 'M', 'baseScan', 'targetScans');
    hi = selectDataType(hi, 'MotionComp_RefScan1');
    saveSession;
    close all;
else 
    warning('Looks between-scan motion was already completed. Will skip this step')
end

fprintf(lid, 'Between-scan motion compensation complete. \n\n');
fprintf('Between-scan motion compensation complete. \n\n');

fprintf(lid, 'Preprocessing for %s is complete! \n', sessionDir);
fprintf('Preprocessing for %s is complete! \n', sessionDir);
fclose(lid);
err = 0;
