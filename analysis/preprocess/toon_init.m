function toon_init(baseDir, expt, session, canonXflag,clip, runBase)
% 
% toon_init(baseDir, expt, session, canonXflag,clip, runBase)
%
% Step 1 of the workflow for analyzing toonotopy
% Initializes the data for mrVista
%
% Default Input values
% baseDir  '/biac2/kgs/projects/'
% expt     'spatiotemporal/experiments/simseq/data'
% session   session name, e.g. 'subj01/session1'
% canonXflag 0; flag to apply cannonical transformation to niftis
% clip      6; number of frames/TRs to clip from beginning of experiment;
% keep      96; number of frames/TRs to keep 
% runBase  'run' ; nifti file tag that designates the toon data
%
% DF 2020 (adapted from code by AS & SP for fLoc)

%% Inputs/Default params
if notDefined('baseDir')
    error('You need to define base directory');
end

if notDefined('session')    
    session = [];
end
if notDefined('canonXflag')
    canonXFlag = 0;
end
if notDefined('clip')
    clip = 5;
end
if notDefined('keep')
    keep = -1; % until the end of TRs (102)
end
if notDefined('runBase')
    runBase = 'run'; % or 'run';
end
if notDefined('rawFileDir')
    rawFileDir = 'raw';
end
%% Checks
% Check and validate inputs and path to vistasoft

if isempty(which('mrVista'))
    vista_path = 'https://github.com/vistalab/vistasoft';
    error(['Add vistasoft to your matlab path: ', vista_path]);
end

% standardize and validate session argument
exptDir = fullfile(baseDir, expt);

dd = dir(exptDir);
allSessions = {dd([dd.isdir]).name};
if length(allSessions) < 3
    error(['No valid session data directories found in ', exptDir]);
    return;
else
    allSessions = allSessions(3:end);
end
if sum(strcmp(session, allSessions)) ~= 1
    error(['Session ', session, ' not found in ', exptDir]);
    return;
end

% look for niftis and create preproc log
sessionDir = fullfile(exptDir, session);
cd(sessionDir);
fns = dir(fullfile(sessionDir,rawFileDir));
fns = {fns.name};
lid = fopen('preprocess_log.txt', 'w+');
fprintf(lid, 'Starting preprocessing analysis for session %s. \n\n', session);
fprintf('Starting preprocessing for session %s. \n\n', session);

%% Optional transformation
% apply cannonical transformation to nifti files and resave
if canonXFlag 
    transformerDir(sessionDir);
end

%% Load niftis and get params
runCount = 0;
for n = 1:length(fns)
    if ~isempty(regexp(fns{n}, runBase))
        if ~contains(fns{n},'.json')
            runCount = runCount + 1;
            niiFiles{runCount} = fullfile(rawFileDir,fns{n});
        end
    end
end
nii = niftiRead(niiFiles{1});
nSlices = size(nii.data, 3);
tr = nii.pixdim(4);
clear nii;

% paths to niftis for each run
niiFiles = natsort(niiFiles);
niiFiles = cellfun(@(X) fullfile(sessionDir, X), niiFiles, 'uni', false);

if length(clip) == 1
    clip = repmat(clip, 1, runCount);
elseif length(clip) ~= length(session)
    fprintf(lid, 'Error -- Length of clip argument is inconsistent with number of runs. \nExited analysis.');
    fprintf('Error -- Length of clip argument is inconsistent with number of runs. \nExited analysis.');
    fclose(lid);
    return;
end
keepFrames = [clip(:), repmat(keep, length(clip), 1)];

%% Initialize mrVista session

% setup analysis parameters 
params = mrInitDefaultParams;
params.doAnalParams = 1;
params.doSkipFrames = 1;
params.doPreprocessing = 0;
params.functionals = niiFiles; % paths to runs of fMRI data
params.subject = session; % name of session directory
params.keepFrames = keepFrames; % TRs to model (after clipping)
params.scanGroups = {1:runCount}; % group all runs of localizer
params.motionComp = 0; % disable motion correction for now
params.sliceTimingCorrection = 1;
params.sliceOrder = [1:2:27,2:2:28]; % Alternating increasing starting at first slice, ending at last

% look for T1 volume and leave blank if none exists
if exist(fullfile(sessionDir, '3DAnatomy', 't1.nii.gz'), 'file') == 2
    params.vAnatomy = fullfile(sessionDir, '3DAnatomy', 't1.nii.gz');
end

% look for inplane volume
dd = dir(fullfile(rawFileDir,'*inplane*'));
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
