function toon_prfRun(baseDir, subjID, session, paramPath, imgPath, prfModel)
% toon_prfRun(baseDir, expt, session);
% 
% Step 4 of the workflow for analyzing toonotopy
% Runs the prf model (CSS)
%
% Example input values
% baseDir       '/oak/stanford/groups/kalanit/biac2/kgs/projects/'
% expt          'subj01'
% session        session name, e.g. 'session1'
% paramPath     'Stimuli/8bars_params.mat';
% imgPath       'Stimuli/8bars_images.mat';
%
% ERK 2021 @ Stanford U, adapted from DF 2020

%% Default inputs
if notDefined('baseDir')
    error('[%s]: You need to define a base directory',mfilename);
end
if notDefined('subjID')
    error('[%s]: You need to define a subject',mfilename);
end
if notDefined('session')
    error('[%s]: You need to define a session',mfilename);
end
if notDefined('paramPath')
    paramPath= 'Stimuli/8bars_params.mat';
end
if notDefined('imgPath')
    imgPath= 'Stimuli/8bars_images.mat';
end
%% Setup
close all; clc;

dataDir = fullfile(baseDir, subjID);
cd(dataDir);

% list of knk scan number, corresponding to session list
scanNum = 1;

% dataTYPE name. Can run for mutiple datatypes
list_rmName = {'Averages'};

% prf model. Specify in a cell. Options:
% {'one oval gaussian' | 'onegaussian' | 'css' | 'difference of gaussians'}
% Note: if we want to specify multiple models, change the naming
% convention. See outFileName
% prfModel = prfModel; %{'onegaussian'};%{'css'};

% search type.
% 1 = grid search only ("coarse"),
% 2 = minimization search only ("fine"),
% 3 = grid followed by minimization search [default]
wSearch = 3;

% radius of circle retinotopy in visual angle degrees
% 20 for wide field and 7 for normal
p.stimSize = 12; % Note S11 has larger stimsize 13.3291;

%% Define things common to all datatypes

% name of params file
p.paramsFile = paramPath;

% image file
p.imFile = imgPath;

% params common to all dts
params.stimSize   = p.stimSize;
params.fliprotate = [0, 0, 0]; %[L/R flip up/down flip rotate degrees]
params.stimType   = 'StimFromScan';
params.stimWidth  = 90; % ??? deg?
params.stimStart  = 0; % Stimulus start volume?
params.stimDir    = 0; % Use alternative stim dir?
params.nCycles    = 1; % Number of stimulus cycles (in case of rotating wedge or expanding ring)
params.nStimOnOff = 0; %??? Number of stimulus on off runs? volumes?
params.nUniqueRep = 1; % number of unique runs
params.nDCT       = 1; % number of DC trends
params.hrfType    = 'two gammas (SPM style)';
params.hrfParams  = {[1.6800, 3, 2.0500], [5.4000, 5.2000, 10.8000, 7.3500, 0.3500]};
params.imfilter   = 'binary';       % convert to binary stimulus
params.jitterFile = 'Stimuli/none'; % ???
params.fixcssexp  = 0;              % Do not fix CSS exponent, but fit it instead
%% Run the prf model
% directory with ret vista session. move here
cd(fullfile(dataDir, session));

% open the session
vw = initHiddenGray;

% need some global variables later
load mrSESSION;

%% loop over the datatypes
for kk = 1:length(list_rmName)

    % set current dataTYPE
    rmName = list_rmName{kk};
    vw = viewSet(vw, 'curdt', rmName);

    % get the dataType struct
    dtstruct = viewGet(vw, 'dtstruct');

    % get the data type number for later
    dataNum = viewGet(vw, 'curdt');

    % ret parameters based on the subject
    % scan number with checkers 
    p.scanNum = scanNum;
    
    % Set parameter and image files
    params.paramsFile = p.paramsFile;
    params.imFile = p.imFile;

    %% getting parameter values for prf model fit ----------------------
    params.nFrames     = viewGet(vw, 'nFrames');
    params.framePeriod = viewGet(vw, 'framePeriod');
    tem.totalFrames    = mrSESSION.functionals(p.scanNum).totalFrames;
    params.prescanDuration = (tem.totalFrames - params.nFrames) * params.framePeriod;

    % store it
    dataTYPES(dataNum).retinotopyModelParams = params;

    % save it
    saveSession;

    %% Put the rm params into the view structure

    vw = rmLoadParameters(vw);
    % the function rmLoadParameters used to call both rmDefineParameters
    % and rmMakeStimulus. If we do it here so that we can give it arguments
    % outside of the default (eg previously, sigma major and minor would be
    % identical despite having prfModel = {'one oval gaussian'} when
    % specifying it as an argument in vw = rmMain(vw, ...)

    % store params in view struct
    vw = viewSet(vw, 'rmParams', params);

    %% RUN THE PRF!

    % name the ret model - whole brain
    outFileName = ['retModel-', prfModel{1}, 'Fit'];

    % no need to load rois, just run it!
    vw = rmMain(vw, [], wSearch, 'model', prfModel, 'matFileName', outFileName);
    
end

