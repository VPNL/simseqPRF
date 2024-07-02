%% s_runDataPreproc.m
% Main script to detrend, average, combine, and chop measured time series.
%
% Written by Eline R Kupers, Stanford U 2021

projectDir = simseqRootPath;
subjnrs    = [1,2,3,7,8,9,10,11,12,13];
roiType0   = 'stimcorner4_area4sq_eccen5';
if useSTRetParams
    temporalModel = '1ch-dcts';  %'1ch-dcts' or '3ch-stLN_opt'; % DN-ST or CSTopt stRet model fit params 
    spatialModel  = 'onegaussianFit'; % referring to stRet  model fit params
    postFix       = []; % '_opt' or '_fix' for CST model
else
    temporalModel = '1ch-glm';
    spatialModel  = 'cssFit';
end

for subjnr = subjnrs
    pths  = getSubjectPaths(projectDir,subjnr);
    
    if useSTRetParams
        roiType = sprintf('%s_stRet_CSTopt_DNST_matchingVoxels',roiType0);
        saveDataFolder = fullfile(pths.simseqResultsDir,...
            sprintf('%s_preprocData_stRet_matchingVoxels', pths.subjID));
        subFolder = sprintf('preprocData_stRet_matchingVoxels_%s%s',temporalModel, postFix);
    else
       roiType = roiType0;
       saveDataFolder = fullfile(pths.simseqResultsDir);
    end
    
    % Define inputs such as folders
    stimFileName   = 'stim_simseq_varyStimDur_x_Area4_2x2design_1deg2_bothf_mid_fullrun';
    removedTRsAtStart = 0; %8;
    simseq_processTimeSeriesData(subjnr, projectDir,...
        'saveDataFolder',saveDataFolder,...
        'subFolder',subFolder, ...
        'roiType',roiType, ...
        'hemi', 'both', ...
        'stimFileName',stimFileName,...
        'bslCorrect', false);
end