%% toon_workflow.m
% Script that runs Toontopy (retinotopy with cartoon images) analysis
%
% Dependencies:
% * vistasoft
% * spm12
% * code folder in Toon experiment folder on native oak:
%   addpath(genpath('/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/toon/code'));
% Add toolboxes with ToolboxToolbox: tbUse({'vistasoft', 'spm12'})
% or manually:
% addpath(genpath('/share/kalanit/software/spm12'));
% addpath(genpath('/share/kalanit/software/vistasoft'));

%% Set session information and paths 
nativeOak           = '/oak/stanford/groups/kalanit/biac2/kgs/';
cd(nativeOak)
nativeOakProjectPth = fullfile(nativeOak,'projects'); 

anatDir = fullfile('./anatomy');
expDir  = fullfile('./projects','spatiotemporal','experiments','toon');
codeDir = fullfile(expDir, 'code');

%% Define subject

subjNumber = 14;
anatNames  = {'kgs201310','kis202008','ek202011', 'hk2019', ... S1-S4
               'lm202105','ac202105','jj202105','jr202105', ... S5-S8
               'xy202105', 'yc202106', 'yf202202','em201611', ... % S9-S12
               'bf202111','st202204'}; % S13, S14,
usedSessions = {'session1','session1','session2', 'session1', ... S1-S4
               'session3','session1','session1','session1', ... S5-S8
               'session1','session1', 'session1','session1_wf', ... % S9-S12
               'session1','session1'}; % S13, S14,
subjID       = sprintf('subj%02d', subjNumber);

anatName    = anatNames{subjNumber};
session     = usedSessions{subjNumber};
sessionPth  = [ 'data/' subjID '/' session];
paramFile   = fullfile('./Stimuli','8bars_params.mat');
imgFile     = fullfile('./Stimuli','8bars_images.mat');


%% Make symbolic link to main anatomy folder
% cd(nativeOakProjectPth)
% vistaPath = fullfile(anatDir,'vistaVol',anatName);gkwlsm
% if ~isfolder(fullfile(expDir, sessionPth, '3DAnatomy'))
%     line = sprintf('ln -s %s %s/3DAnatomy',vistaPath,fullfile(expDir, sessionPth));
%     system(line);
% end

cd(nativeOakProjectPth)
if ~exist(fullfile(expDir, sessionPth, '3DAnatomy'),'dir')
    vistaPath = fullfile('../../../../../../../',anatDir,'vistaVol',anatName);
    line = sprintf('ln -s %s %s/3DAnatomy',vistaPath,fullfile(expDir, sessionPth));
    system(line);
end
%% Set up anatomy path
cd(fullfile(nativeOakProjectPth,expDir,sessionPth)); %fix the path
anatSubjPth = fullfile('./3DAnatomy','t1.nii.gz');
setVAnatomyPath(anatSubjPth);
% saveSession;

%% Copy stimfile and remove color stimulus information (ONLY NEEDED ONCE)
yFlip = 0;
if exist('./Stimuli/8bars_images.mat', 'file') ~= 2
    copyfile ../../subj01/session1/Stimuli/ ./Stimuli/
    load('./Stimuli/8bars_images.mat')
    if length(size(images)) ~=3
        for i = 1:size(images,4)
            bk_img(:,:,i) = rgb2gray(images(:,:,:,i));
        end
        images = bk_img;
        save('./Stimuli/bars_9steps.mat','images');
        clear bk_img;
    end
    
    if yFlip ==1 % if you need to flip y axis
        for oneFrame = 1:size(images,3)
            images(:,:,oneFrame) = flipud(images(:,:,oneFrame));
            images(:,:,oneFrame) = fliplr(images(:,:,oneFrame));
        end
        save('./Stimuli/8bars_images.mat','images');
    end
end
%% Initialize session (ONLY NEEDED ONCE)
if ~exist(fullfile(nativeOakProjectPth,expDir, ['data/' subjID], session,'Inplane'),'dir')
    toon_init(fullfile(nativeOakProjectPth,expDir),  ['data/' subjID], session);
end

% To check initialization
% cd to subject's directory to initialize mrVista (type mrVista)
% Using GUI: load mean map & check that functionals match Inplane anatomicals
% Best visualization: threshold mean map at around 500-1000 to see anatomical underlay

%% Align inplane anatomy to volume anatomy (ONLY NEEDED ONCE)
load('mrSESSION.mat')
if ~isfield(mrSESSION, 'alignment') || isempty(mrSESSION.alignment)
%     rxAlign;
%     s_alignInplaneToAnatomical
end

%% motion correct session
if ~exist(fullfile('./Inplane/MotionComp_RefScan1'),'dir')
    toon_motionCorrect(fullfile(nativeOakProjectPth,expDir, 'data'), subjID, session);
end
%% install segmentation, transform tSeries to Gray, and average time series
% Use getVAnatomyPath to check VAnatomyPath
setVAnatomyPath(anatSubjPth)
if ~exist(fullfile('./Gray/MotionComp_RefScan1'),'dir')
    toon_2gray(fullfile(nativeOakProjectPth,expDir,'data'), subjID, session)
end

%% run CSS pRF model (to run other models, see line 48 in toon_pRFRun.m)
% Specify as cell, for example: {'onegaussian'}.
% You can run multiple models at once, for example: {'onegaussian', 'css'}.
% Choose from: 
%       'onegaussian'      :    your standard vanilla 2D circular Gaussian
%       'one oval gaussian':    eliptical pRF
%       'css'              :    compressive spatial summation pRF 
%                               (see Kay et al. 2013 J.Neurophys)
%       'difference of gaussians': Having an positive 2D Gaussian and
%                               equally wide or wider, negative 2D Gaussian
%                               (see Zuiderbaan et al.2012 JoV)
prfModel = {'css'}; %{'css'};
toon_prfRun(fullfile(nativeOakProjectPth,expDir,'data'), subjID, session, paramFile, imgFile, prfModel)
  