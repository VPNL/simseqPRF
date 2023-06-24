function simseq_init(sessionDir,clip, keep, runBase, rawFileDir)
% simseq_init(sessionDir, clip, keep, runBase, rawFileDir)
%
% Function to initialize workflow for mrVista
%
% Inputs:
% sessionDir        : (str) base folder of session, where mrSESSION.mat
%                       will be created
% clip              : (int) nr of frames/TRs clipped from every functional run
% keep              : (int) number of frames/TRs to keep, use -1 for all.
% runBase           : (str) common string in EPI file name to identify files
% rawFileDir        : (str) folder where raw data live
%
% Written by EK 2021/01 Stanford U (adapted from toonotopy code by D Finzi)

%% Find session path and raw data

[filepath, ~] = fileparts(sessionDir);
[~, subjID]   = fileparts(filepath);

cd(sessionDir);
fns = dir(fullfile(sessionDir,rawFileDir));
fns = {fns.name};

%% Start preprocess log
lid = fopen('preprocess_log.txt', 'w+');
fprintf(lid, 'Starting preprocessing analysis for session %s. \n\n', sessionDir);
fprintf('Starting preprocessing for session %s. \n\n', sessionDir);

%% Load niftis and get params
runCount = 0;
for n = 1:length(fns)
    if contains(fns{n}, runBase)
        runCount = runCount + 1;
        niiFiles{runCount} = fullfile(rawFileDir,fns{n});
    end
end

% Get params
nii     = niftiRead(niiFiles{1});
nSlices = size(nii.data, 3);
tr      = nii.pixdim(4);
clear nii;

% Paths to niftis for each run
niiFiles = natsort(niiFiles);
niiFiles = cellfun(@(X) fullfile(sessionDir, X), niiFiles, 'uni', false);

if length(clip) == 1
    clip = repmat(clip, 1, runCount);
elseif length(clip) ~= runCount
    fprintf(lid, 'Error -- Length of clip argument is inconsistent with number of runs. \nExited analysis.');
    fprintf('Error -- Length of clip argument is inconsistent with number of runs. \nExited analysis.');
    fclose(lid);
    return;
end
keepFrames = [clip(:), repmat(keep, length(clip), 1)];

%% Initialize mrVista session

% setup analysis parameters 
params = mrInitDefaultParams;
params.doAnalParams    = 1;
params.doSkipFrames    = 1;
params.doPreprocessing = 0;
params.functionals     = niiFiles;      % paths to runs of fMRI data
params.subject         = subjID;        % name of session directory
params.keepFrames      = keepFrames;    % TRs to model (after clipping)
params.scanGroups      = {1:runCount};  % group all runs of localizer
params.motionComp      = 0;             % disable motion correction for now
params.sliceTimingCorrection = 1;
params.sliceOrder = [1:2:16,2:2:16]; % Alternating increasing starting at first slice, ending at last

% Look for T1 volume and leave blank if none exists
if exist(fullfile(sessionDir, '3DAnatomy', 't1.nii.gz'), 'file') == 2
    params.vAnatomy = fullfile(sessionDir, '3DAnatomy', 't1.nii.gz');
end

% look for inplane volume
dd = dir('*Inplane*');
if isempty(dd) 
    dd = dir(fullfile(rawFileDir,'*Inplane*'));
end

inplane = dd.name;
if isempty(inplane)
    fprintf(lid, 'Warning -- Inplane scan not found. Continued analysis. \n');
    fprintf('Warning -- Inplane scan not found. Continued analysis. \n');
else
    params.inplane = fullfile(sessionDir,rawFileDir,inplane);
end

% inititalize vistasoft session and open hidden inplane view
fprintf(lid, 'Initializing vistasoft session directory in: \n%s \n\n', sessionDir);
fprintf('Initializing vistasoft session directory in: \n%s \n\n', sessionDir);
if ~exist(fullfile(sessionDir, 'Inplane'), 'dir')
    mrInit(params);
end
end