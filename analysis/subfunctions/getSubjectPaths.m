function pths = getSubjectPaths(projectDir, subjID)
% Function to get subject path to reproduce data/figures of paper:
%
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University

%%
pths = struct();
pths.projectDir       = projectDir;
pths.subjnr           = subjID;
pths.subjID           = sprintf('subj%02d', subjID);
% pths.anatName       = allAnatNames{subjID}; % no anatomy in this folder
pths.dataDirSimSeq    = fullfile(projectDir, 'data', 'simseq');
pths.dataDirToon      = fullfile(projectDir, 'data', 'toon');
pths.stimDir          = fullfile(projectDir, 'data', 'stimuli');
pths.figureDir        = fullfile(projectDir, 'results',pths.subjID,'figs');
pths.simseqResultsDir = fullfile(projectDir, 'results',pths.subjID);

if subjID == 12
    pths.toon   = fullfile(pths.dataDirToon,sprintf('%s/%s',pths.subjID,'session1_wf'));
else
    pths.toon   = fullfile(pths.dataDirToon,sprintf('%s/%s',pths.subjID,'session1'));
end

% SIMSEQ SUBJECT SESSIONS
switch subjID
    case 1 % subj01
        pths.expversionNr  = 11; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 2 % subj02
        pths.expversionNr  = 5; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 3 % subj03
        pths.expversionNr  = 10; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 7 % subj07
        pths.expversionNr  = 7; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 8 % subj08
        pths.expversionNr  = 14; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 9 % subj09
        pths.expversionNr  = 6; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 10 % subj10
        pths.expversionNr  = 8; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2];
        pths.scansToLoad   = {[1:6]};
    case 11 % subj11
        pths.expversionNr  = 9; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 12 % subj12
        pths.expversionNr  = 12; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
    case 13 % subj13
        pths.expversionNr  = 15; % has 333 TRs
        pths.runOrder      = [1 2 1 2 1 2 1 2];
        pths.scansToLoad   = {[1:8]};
end

% %% Stimuli
d = dir(fullfile(pths.dataDirSimSeq,'behavior',pths.subjID,'simseqPRF*.mat'));

for ii = 1:length(d)
    pths.stimFiles{ii} = fullfile(d(ii).folder,d(ii).name);
end

%% ROIs
% We exclude ROIs for some subjects because they are not retinotopically
% defined (for example, voxels falls outside of the FoV of the EPI)
pths.preferredRoiName = @(h,r) sprintf('%s_%s_toon',h,r);
pths.allROIs = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', ...
    'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO','LO','TO'};

switch subjID
    case 1
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2','TO1', 'TO2', 'IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2', 'TO1', 'IPS0'};
        pths.definedROIsBOTH = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2', 'TO1', 'IPS0'};
    case 2
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2','TO1', 'TO2','IPS0','IPS1'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2','V3AB','LO1','LO2','TO1', 'TO2'};
        pths.definedROIsBOTH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2','V3AB','LO1','LO2','TO1'};
    case 3
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB', 'LO1', 'LO2', 'TO1','TO2','IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB', 'LO1', 'LO2', 'TO1'};
        pths.definedROIsBOTH = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB', 'LO1', 'LO2', 'TO1'};
    case 7
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2'};
        pths.definedROIsBOTH = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2'};
    case 8
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0','IPS1'};
        pths.definedROIsBOTH  = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0','IPS1'};
    case 9
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2','IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','IPS0','IPS1'};
        pths.definedROIsBOTH = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','IPS0','IPS1'};
    case 10
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','TO1','TO2', 'IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1', 'V3AB','LO1','LO2','TO1','TO2'};
        pths.definedROIsBOTH = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2','TO1','TO2'};
    case 11
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2','TO1'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1'};
        pths.definedROIsBOTH  = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1'};
    case 12
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO2','TO1','TO2'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2'};
        pths.definedROIsBOTH  = {'V1', 'V2', 'V3', 'hV4', 'VO1','V3AB','LO1','LO2'};
    case 13
        pths.definedROIsLH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0'};
        pths.definedROIsRH   = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0'};
        pths.definedROIsBOTH  = {'V1', 'V2', 'V3', 'hV4', 'VO1','VO2', 'V3AB','LO1','LO2','TO1','TO2', 'IPS0'};
end


return

