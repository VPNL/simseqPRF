%% s_runDataPreproc.m
% Main script to detrend, average, combine, and chop measured time series.
%
% Subjnrs: [1,2,3,7,8,9,10,11,12,13]
% Written by Eline R Kupers, Stanford U 2021

projectDir = simseqRootPath;

for subjnr = 3
    sesNr = getSessionNrMainExp(subjnr);
    pths  = getSubjectPaths(projectDir,subjnr,sesNr);
    
    % Define inputs such as folders
    stimFileName   = 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun';
    subFolder      = 'varyStimDur_x_Area4_2x2design_4deg2_trimmedGaussian_bothHemi_mid_roiCorrected_fullrun';
    saveDataFolder = fullfile(pths.simseqResultsDir,'preprocData');
    removedTRsAtStart = 0; %8;
    simseq_processTimeSeriesData(subjnr, projectDir,...
        'saveDataFolder',saveDataFolder,'subFolder',subFolder, 'sessionNr', sesNr, ...
        'roiType','stimcorner4_area4sq_eccen5', 'hemi', 'both', 'stimFileName',stimFileName,...
        'bslCorrect', false);
end