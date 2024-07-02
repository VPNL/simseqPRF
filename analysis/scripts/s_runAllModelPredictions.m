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
roiType             = 'stimcorner4_area4sq_eccen5';
spatialModels       = {'onegaussianFit','cssFit','onegaussianFit'}; % choose from: 'onegaussianFit','cssFit'
temporalModels      = {'1ch-glm','1ch-glm','3ch-stLN'}; % choose from: '1ch-glm','3ch-stLN','1ch-dcts'
useSTRetParams      = true; % If true, we use solved parameters from spatiotemporal retinotopy paper for 3ch-stLN model

if useSTRetParams
    subjnrs = [1:3,7:10];
    subDir = 'savedPredictions_fixedPRF_combT_matchingVoxels';
else
    subjnrs = [1:3,7:13];
    subDir = 'savedPredictions_fixedPRF_combT';
end

if verbose, makeprettyfigures; end %#ok<UNRCH>

for subjnr = subjnrs
    
    pths = getSubjectPaths(projectDir, subjnr);
    simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);
    versionNr        = pths.expversionNr;
    
    for runNr = runNrs
        
        % --------------------------------------
        % -- Get experiment folders and files --
        % --------------------------------------
        saveFolder = sprintf('%s_fullrun%d_v%d', runNr, versionNr);
        savePredFolder0 = fullfile(pths.simseqResultsDir,subDir,saveFolder);
        fname = sprintf('stim_simseq_run%d_v%d.mat',runNr, versionNr);
        stimfile = fullfile(simSeqSessionDir, 'Stimuli', fname);
        % -------------------------------------
        
        for ii = 1:length(spatialModels)
            
            if strcmp(temporalModels{ii},'1ch-glm') && ~useSTRetParams
                % Standard Gaussian model
                useFixedExponent = [];
                exponents        = NaN;
                combChan         = [];
                savePredFolder   = fullfile(savePredFolder0);
                roiType          = roiType0;
                
            elseif strcmp(temporalModels{ii},'1ch-dcts') && useSTRetParams
                % Standard Gaussian model + 1-chan Divisive Normalization spatiotemporal model
                useFixedExponent = [];
                combChan         = [];
                exponents        = NaN;
                savePredFolder   = fullfile(savePredFolder0,sprintf('modelPredictions_stRet_%s',temporalModels{ii}));
                roiType          = sprintf('%s_%s',roiType0,'CSTopt_DNST_matchingVoxels');
                
            elseif strcmp(temporalModels{ii},'3ch-stLN') && ~useSTRetParams
                % Standard Gaussian model + 3-chan linear-nonlinear spatiotemporal model
                combChan         = [1 2 2];
                savePredFolder   = fullfile(savePredFolder0,'modelPredictions_gridfit2');
                exponents        = [0.1:0.05:1]; % choose any value between 0-1 or NaN for no fixed exponent in modelfit
                roiType          = roiType0;
                
            elseif strcmp(temporalModels{ii},'3ch-stLN') && useSTRetParams
                % Standard Gaussian model + 3-chan linear-nonlinear
                % spatiotemporal model with parameters from ST ret exp
                combChan         = [1 2 2];
                savePredFolder   = fullfile(savePredFolder0,sprintf('modelPredictions_stRet_%s',temporalModels{ii}));
                exponents        = NaN; % We'll use the parameters from the spatiotemporal retinotopy experiment (Kim et al. 2024)
                roiType          = sprintf('%s_%s',roiType0,'CSTopt_DNST_matchingVoxels');
            end
            
            
            for expn = exponents
                
                simseq_predictTimeseriesFromToonPRFs(subjnr, projectDir, stimfile,...
                    'spatialModel',spatialModels{ii},...
                    'temporalModel',temporalModels{ii}, ...
                    'upsampleStimType',upsampleStimType,...
                    'saveFolder',savePredFolder, ...
                    'veThresh',veThresh, ...
                    'verbose', verbose, ...
                    'hemi',hemi, ...
                    'roiType',roiType, ...
                    'stimRun',runNr, ...
                    'useArtificialPRFs',false, ...
                    'useMedianROIexponent', false, ...
                    'combineNeuralChan', combChan, ...
                    'useFixedExponent',expn, ...
                    'useSTRetParams',useSTRetParams,...
                    'roiIdx','all');
            end
        end
    end
end
    
return
