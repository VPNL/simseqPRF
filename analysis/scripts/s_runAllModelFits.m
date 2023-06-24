%% s_runAllModelFits.m
% Main script to fit predict time series to measured time series.
%
% Some notes about this script:
% * previous version of this function is:
% simseq_fitTemporalModelPredictiosToDataWithPreproc('stimFileName',stimFileName, ...)
% * We remove 10 TRs at the start of each predicted run's time series and
% the stimulus sequence, (not the preprocessed data, see s_runDataPreproc.m)
% because the data preprocessing removed 6 TRs (1 TR = 1 sec) countdown
% + 2 TRs of the baseline. Here we remove 10TRs, so two extra TRs, because
% we observed that there is still gradient instability at the start of each
% time series run. Keeping this introduces artefacts in baseline correction
% and results in poorer modelfits.
%
% Main experiment subject numbers = [1,2,3,7,8,9,10,11,12,13]
%
%
% Written by Eline R Kupers, Stanford U 2021

%% Define paths and parameters
projectDir = '/oak/stanford/groups/kalanit/biac2/kgs/projects/spatiotemporal/';
pilotNr    = 3;
exponents  = [NaN,NaN, 0.1:0.05:1]; %Use NaN for linear and CSS model, use [0.1:0.1:1] for ST grid fit

for subjnr = [1,2,3,7,8,9,10,11,12,13]
    
    sesNr = getSessionNrMainExp(subjnr);
    pths  = getSubjectPaths(projectDir,subjnr,sesNr);
    
    % Define inputs such as folders
    stimFileName   = 'stim_simseq_run';
    subFolder      = 'varyStimDur_x_Area4_2x2design_fullrun';
    saveBetaFolder = fullfile(pths.simseqResultsDir,'savedPredictions_fixedPRF_combT');
    loadDataFolder = fullfile(pths.simseqResultsDir,'preprocData');
    saveFigDir     = fullfile(pths.figureDir, 'cvfits');
    removedTRsAtStart = 2; % Preprocessing EPIs removed 8 TRs (= 6 TRs countdown + 2 TRs baseline),
    % Here, we remove two extra TRs, since the countdown is not part of model prediction
    spatialModels   = {'onegaussianFit','cssFit','onegaussianFit'};
    temporalModels  = {'1ch-glm','1ch-glm','3ch-stLN'};
    
    for expn = exponents
        if isnan(expn)
            modelsToRun = [1,2];
            useFixedExponent = NaN;
        else
            modelsToRun = 3;
            useFixedExponent = expn;
        end
        for ii = modelsToRun
            simseq_fitTemporalModelPredictionsToPreprocData(...
                subjnr, projectDir,spatialModels{ii},temporalModels{ii},...
                'saveBetaFolder',saveBetaFolder,'subFolder',subFolder, ...
                'sessionNr', sesNr,'saveFigDir',saveBetaFolder, ...
                'loadDataFolder',loadDataFolder, ...
                'roiType','stimcorner4_area4sq_eccen5', ...
                'hemi', 'both', ...
                'removedTRsAtStart',removedTRsAtStart, ...
                'weightChannels', [], ...
                'doCrossValidation', true, ...
                'roiIdx','all', ...
                'useFixedExponent',useFixedExponent,...
                'useSearchFit',false, ...
                'verbose', true);
        end
    end
    
end

end