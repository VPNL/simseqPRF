%% simseq_preprocess_workflow_raw_to_gray.m
%
% Script to preprocess sim seq pilot data.
% Includes the following steps:
% - Initialize (w/ slice timing correction)
% - Motion correct
% - Align inplane to T1w anatomy
% - Convert voxel data to mrVista gray nodes and average/group data
% 
% Make sure to remove symbolic toolbox from path, as it overwrites
% vistasoft functions such as unitConvert
% rmpath(genpath('/software/matlab/R2020b/toolbox/symbolic/'))
%
% Toolbox dependencies:
% * knkutils
% * alignvolumedata
% * vistasoft
% * simseqPRF
%
% use tbUse({'simseqPRF','vistasoft','alignvolumedata','knkutils'})
% or
% addpath(genpath('~/matlab/git/toolboxes/simseqPRF/'))
% addpath(genpath('~/matlab/git/toolboxes/knkutils/'))
% addpath(genpath('~/matlab/git/toolboxes/alignvolumedata/'))
% addpath(genpath('/share/kalanit/software/vistasoft/'))

% Written by ERK @ Stanford U 2021

%% Define subject

% cd to subject's session folder.
cd /oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/experiments/simseq/

projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal';
anatDir    = '/oak/stanford/groups/kalanit/biac2/kgs/anatomy';
subjnr     = 11;
sessionNr  = getSessionNrMainExp(subjnr);
pths       = getSubjectPaths(projectDir, subjnr,sessionNr);

cd(fullfile(pwd, 'data', pths.subjID, pths.session));

%% Initialize vistasoft session 
% needed for motion correction and inplane alignment

sessionDir = pwd;
clip       = 8;    % 8 seconds = 8 TRs of beginning  (previous pilot 1 and 2 had 6 TRs removed)
keep       = -1;   % use all frames until the end
runBase    = 'run';
rawFileDir = 'raw';
simseq_init(sessionDir,clip, keep, runBase, rawFileDir);

%% Within and between scan motion correction
simseq_motionCorrect(sessionDir);

%% Set up 3D anatomy path
if ~isfolder(fullfile(pths.dataDirSimSeq, pths.subjID,pths.session,'3DAnatomy'))
    vistaPath = fullfile(anatDir,'vistaVol',pths.anatName);
    line = sprintf('ln -s %s %s/3DAnatomy',vistaPath,...
        fullfile(pths.dataDirSimSeq, pths.subjID,pths.session));
    system(line);
end

anatSubjPth = fullfile('./3DAnatomy','t1.nii.gz');
setVAnatomyPath(anatSubjPth);

%% Align inplane to T1 anatomy with rxAlign
%% NB: use  addpath(genpath('/software/matlab/vistasoft/'))
% This step has to happen after simseq_init and before simseq_2gray,
% but can happen before or after simseq_motionCorrect.
% Use:
% s_alignInplaneToAnatomical 
%
% NB. this requires knkutils and alignvolumedata toolboxes by Kendrick Kay,
% Script starts with rxAlign to get inplane roughly in the right place,
% rest should be pretty automatic (except for defining elipse to set a
% roughly estimated brain/skull boundary)

%% Convert voxel data to gray nodes
simseq_2gray(sessionDir)

%% Add par files to session
% Here we do it by script (see er_assignParfilesToScans), which is similar
% to doing it in the GUI: mrVista 3
%
% Select MotionComp_RefScan1
% Go to GLM > assign parfiles > select run 1:2:end and assign to run 1
% Go to GLM > assign parfiles > select run 2:2:end and assign to run 2
% Go to GLM > grouping > group all scans > select MotionComp_RefScan1 >
% select all scans.
% Go to File > Save mrSESSION
% load('mrSESSION.mat')
% hg = initHiddenGray('MotionComp_RefScan1', 1);
for dt = 1:length(dataTYPES)
    if strcmp(dataTYPES(dt).name,'MotionComp_RefScan1')
        curdt = dt;
        break;
    end
end

nrScans = length(dataTYPES(curdt).scanParams);
d = dir(fullfile('./Stimuli/parfiles/simseq*_run1_*.par'));
for scan = 1:2:nrScans
    dataTYPES(curdt).scanParams(scan).parfile = fullfile(d.name);
end

d = dir(fullfile('./Stimuli/parfiles/simseq*_run2_*.par'));
for scan = 2:2:nrScans
    dataTYPES(curdt).scanParams(scan).parfile = fullfile(d.name);
end
save('mrSESSION', 'dataTYPES', '-append');

% group scans
grpTxt = sprintf('%s: %s',dataTYPES(curdt).name,num2str(1:nrScans));
for s = 1:nrScans
    dataTYPES(curdt).scanParams(s).scanGroup = grpTxt;
end
save('mrSESSION', 'dataTYPES', 'mrSESSION','vANATOMYPATH');

%% Prepare stimulus for event related analysis (only for pilot 1)
% simseq_prepareStimulus(sessionDir)



