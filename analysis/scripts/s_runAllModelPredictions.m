%% s_runAllModelPredictions
% Main script to predict time series for entire run, based on pRF
% parameters from a separate retinotopy experiment. 
% 
% Written by Eline R Kupers, Stanford U 2021

projectDir = simseqRootPath;
verbose    = false; % Print figures with predictions or not.
veThresh   = 0.1; % Variance explained threshold for pRF model estimates
exponents  = [0.1:0.05:1]; % choose any value between 0-1 or NaN for no fixed exponent in modelfit 
runNrs     = [1,2]; 


% -------------- Main experiment params -----------------------
trialType           = 'varyStimDur_x_Area4_2x2design'; 
stimAreaDeg2        = 4; % deg2
hemi                = 'both';
upsampleStimType    = 'upsampleFrom60Hz'; % choose from: 'upsampleFrom60Hz', 'upsampleFrom1Hz', or 'none'
roiType             = 'stimcorner4_area4sq_eccen5'; % choose from: stimcorner4_area2sq_eccen5, stimcorner4_area4sq_eccen5
spatialModels       = {'onegaussianFit','cssFit','onegaussianFit'};
temporalModels      = {'1ch-glm','1ch-glm','3ch-stLN'};

if verbose, makeprettyfigures; end %#ok<UNRCH>

for subjnr = 13
    
    pths = getSubjectPaths(projectDir, subjnr);
    
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session); 
    versionNr        = pths.expversionNr;

    for runNr = runNrs
        
        % --------------------------------------
        % -- Get experiment folders and files --
        % --------------------------------------
        saveFolder = sprintf('%s_fullrun%d_v%d', runNr, versionNr);
        savePredFolder0 = fullfile(pths.simseqResultsDir,'savedPredictions_fixedPRF_combT',saveFolder);
        fname = sprintf('stim_simseq_run%d_v%d.mat',runNr, versionNr);
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', fname);
        % -------------------------------------

        
        for expn = exponents
            if isnan(expn)
                modelsToRun = [1,2];
                useFixedExponent = [];
                savePredFolder = fullfile(savePredFolder0);
            else
                % Standard Gaussian model + 3-chan linear-nonlinear spatiotemporal model
                modelsToRun = 3;
                useFixedExponent = expn;
                savePredFolder = fullfile(savePredFolder0,'modelPredictions_gridfit');
            end
            for ii = modelsToRun
                if ii == 3, combChan = [1 2 2]; 
                else combChan = []; end
                simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, stimfile,...
                    'spatialModel',spatialModels{ii},'temporalModel',temporalModels{ii}, ...
                    'upsampleStimType',upsampleStimType,'saveFolder',savePredFolder, ...
                    'veThresh',veThresh, 'verbose', verbose, ...
                    'hemi',hemi, 'roiType',roiType, 'stimRun',runNr, ...
                    'useArtificialPRFs',false, ...
                    'useMedianROIexponent', false, ...
                    'combineNeuralChan', combChan, ...
                    'useFixedExponent',expn, ...
                    'roiIdx','all');
            end
        end
    end
end

return
