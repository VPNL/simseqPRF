% s_simulateDoGPRFTimeSeries.m

% Set params
projectDir    = fullfile(simseqRootPath);
subjnrs       = [1,2,3,7:13];
hemi          = {'both'};
prfModel      = 3; % Use linear pRF model to get center size
veThresh      = 0.1;
plotFigures   = true;
loadStoredSimulation = false;
roiType       = 'stimcorner4_area4sq_eccen5';
conditionNamesSimSeq = {'Small & short (0.2s)','Small & long (1s)',...
                        'Big & Short (0.2s)','Big & Long (1s)'};

% ROI size ratios from Aqil et al. 2021 PNAS
rois = {'V1','V2','V3','hV4'};
sizeRatio     = [7.4, 6.8, 7.3, 5.8]; % V1,V2,V3,hV4
surroundFun   = @(sigma, roiIdx) sizeRatio(roiIdx).*sigma;

% Load gridfit data
gridFit = load(fullfile(fullfile(projectDir,'experiments/simseq/results/average/gridFit2_Results',...
    'S1_S2_S3_S7_S8_S9_S10_S11_S12_S13_BlockAmps_pRFparams_modelFit_results_cv_variableBlockOnset_gridFit2_v4.mat')));

% Loop over subjects
for subjnr = subjnrs
    s = (subjnr==subjnrs);
    
    % Reset spatial and temporal pRF model
    params.analysis.spatialModel = 'differenceOfGaussiansFit';
    params.analysis.temporalModel = '1ch-glm';
    
    % Define subject paths
    pths   = getSubjectPaths(projectDir, subjnr);
    simSeqDataDir = fullfile(pths.dataDirSimSeq, pths.subjID);
    cd(simSeqDataDir)
    saveFolder = fullfile(pths.simseqResultsDir,'DogSimulation');
    saveFigureDir = fullfile(pths.figureDir,'DogSimulation');
    if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
    
    % Define analysis params
    params.analysis.temporal.param.exponent = [];
    params.analysis.temporal.param.fs = 1000; % hz
    params.analysis.temporal.fs = 1000; % hz
    params.analysis.temporal.tr = 1; % s
    params.analysis.normAcrossRuns = false;
    params.analysis.hrf.func = [];
    params.analysis.hrf.type = 'spm';
    params.useGPU = false;
    params.saveDataFlag = true;
    params.analysis.normNeuralChan        = false;
    params.stim.framePeriod               = 1;
    params.saveDataFlag                   = true;
    params.analysis.spatial.sparsifyFlag  = false; % leave out zero response in conv(prf,stim)
    params.analysis.spatial.keepPixels    = [];
    params.recomputePredictionsFlag       = true;
    params.analysis.zeroPadPredNeuralFlag = true;
    params.analysis.reluFlag              = false;
    params.analysis.useMedianROIexponent  = false;
    
    % Add XY fov grid
    params.analysis.stim.nSamples   = 50; % 50 corresponds to 101 x 101
    params.analysis.stim.fieldSize  = 12; % deg visual angle
    params.analysis.stim.sampleRate = params.analysis.stim.fieldSize/params.analysis.stim.nSamples; % deg visual angle
    x = single( linspace(- params.analysis.stim.fieldSize,  params.analysis.stim.fieldSize, 1+ params.analysis.stim.nSamples*2) );
    y = x;
    [X,Y] = meshgrid(x,y);
    
    % Update the sampling grid to reflect the sample points used.
    params.analysis.spatial.X = X(:);
    params.analysis.spatial.Y = -1*Y(:);
    params.analysis.numberStimulusGridPoints = params.analysis.stim.nSamples;
    params.analysis.sampleRate = params.analysis.stim.sampleRate;
    params.analysis.fieldSize  = params.analysis.stim.fieldSize;
    
    %%
    for roiIdx = 1:length(rois)
        
        if loadStoredSimulation
            load(fullfile(simSeqDataDir, 'DogSimulation', ...
                sprintf('dogSimulation_%s_%s_%s.mat',pths.subjID, rois{roiIdx}, roiType)));
            stim = logical(params.stim.images);
            
            loadFile = fullfile(simSeqDataDir,...
                sprintf('preprocData_%s_%s_%s_%s.mat', ...
                rois{roiIdx}, roiType,'cssFit', 'run12'));
            load(loadFile, 'T','allDataTrials','data','trials',...
                'nrTrials','er_data','stimMS','alignMe');
            data.predictionCatRuns = cat(1,stPred{1}.predBOLD(3:end,:),stPred{2}.predBOLD(3:end,:));
            
        else
            loadFile = fullfile(simSeqDataDir, ...
                sprintf('%s_%s_onegaussianFit_3ch-stLN_cvfitResults', ...
                pths.subjID, rois{roiIdx}));
            gridFit = load(loadFile);
            
            params.analysis.spatial.x0 = gridFit.data.params.x0.both;
            params.analysis.spatial.y0 = gridFit.data.params.y0.both;
            params.analysis.spatial.effectiveSize = gridFit.data.params.effectiveSize.both;
            params.analysis.spatial.exponent =  gridFit.data.params.exponent.both;
            params.analysis.spatial.exponent_temporal = gridFit.data.params.exp_temporal.both;
            params.analysis.spatial.sigmaMajor = params.analysis.spatial.effectiveSize;
            params.analysis.spatial.sigmaMinor = params.analysis.spatial.effectiveSize;
            params.analysis.spatial.sigmaSurround = surroundFun(params.analysis.spatial.effectiveSize,roiIdx);
            
            stPred = cell(1,2);
            for rr = 1:2
                
                % Get stimulus file
                fname    = sprintf('stim_simseq_run%d_v%d.mat', rr, pths.expversionNr);
                stimfile = fullfile(pths.stimDir,pths.subjID, fname);

                params.analysis.predFile = fullfile(saveFolder,  ...
                    sprintf('modelPredictions_%s_%s_%s_run%d.mat',...
                    roi,roiType, params.analysis.spatialModel, params.analysis.temporalModel,rr));
                
                % Go to subjects session directory
                cd(simSeqSessionDir)
                
                % Add stim params
                params.stim.framePeriod     = 1; %s (sample rate of Simseq fMRI data)
                params.stim.imFile          = stimfile;
                params.stim.paramsFile      = stimfile;
                params.stim.prescanDuration = 0; % prescan duration is already clipped during preprocessing
                params.analysis.pRFmodel = {'st-upsampling60hz'};
                
                % Make 1ms stimulus
                params = makeStiminMS(params);
                stim   = logical(params.stim.images);
                tmp    = stPredictBOLDFromStim(params, stim);
                stPred{rr} = tmp;
            end
            
            loadFile = fullfile(simSeqDataDir,...
                sprintf('preprocData_%s_%s_%s_%s.mat', ...
                rois{roiIdx},roiType,'cssFit', 'run12'));
            load(loadFile, 'T','allDataTrials','data','trials',...
                'nrTrials','er_data','stimMS','alignMe');
           
            er_data(1).params.timeWindow = repmat({-4:1:18}, length(trials.condNums), 1);
            er_data(1).params.bslPeriod  = repmat({-4:1:0}, length(trials.condNums), 1);
            er_data(1).params.normBsl    = 0;
            er_data(1).params.parfiles   = trials.parfiles;
            
            data.predictionCatRuns = cat(1,stPred{1}.predBOLD(3:end,:),stPred{2}.predBOLD(3:end,:));
            
            params.analysis.regressionType = 'OLS';
            Y = data.meanDetrendedData;
            X = data.predictionCatRuns;
            [nFrames, numVoxels] = size(data.meanDetrendedData);
            
            % Split data in even and odd concatenated runs
            Y_crossval = {data.meanDetrendedDataOdd, data.meanDetrendedDataEven};
            
            % If we use fracridge, use best alpha from full fit.
            bestAlpha = [];
            
            % Loop over split halves
            for cx = [1,2]
                resultsCValInitFit{cx} = fitModelToVoxData(...
                    X,Y_crossval{cx},...
                    'regressionType',params.analysis.regressionType, ...
                    'sumSTFlag', false, ...
                    'weightChannels', [], ...
                    'addOffsetFlag', true, ...
                    'alpha',bestAlpha);
            end
            
            % Cross-validate modelfits, i.e. use scaled modelfit from split A
            % calculate variance explain with data in split B and vice versa)
            B_crossval = NaN(1,numVoxels);
            R2_crossval = NaN(2, numVoxels);
            for cx = [1,2]
                % We refit beta weights to data where baseline is removed
                resultsCValFinalFit{cx} = fitModelToVoxData(...
                    X,resultsCValInitFit{cx}.Y_noOffset,...
                    'regressionType',params.analysis.regressionType, ...
                    'sumSTFlag', false, ...
                    'weightChannels', [], ...
                    'addOffsetFlag', false, ...
                    'alpha',bestAlpha);
                B_crossval(:,:,cx) = resultsCValFinalFit{cx}.B;
                
                % Get opposite dataset
                fitIdx = setdiff([1,2],cx);
                
                % First, get prediction and data final fit, no offset subtracted
                Xcx = resultsCValFinalFit{cx}.sumChannelPrediction;
                Ycx = resultsCValInitFit{fitIdx}.Y_noOffset;
                
                % Now compute Coefficient of Determination (R2) using split
                % half data B, and modelfit A, and vice versa.
                for n = 1:size(Ycx,2)
                    R2_crossval(cx,n) = computeCoD(Ycx(:,n),Xcx(:,n));
                end
                
                % Also store predictions of final fit with, and without offset
                sumChannelPredictionCVNoOffset{cx} = resultsCValFinalFit{cx}.sumChannelPredictionNoOffset;
                sumChannelPredictionCV{cx}         = resultsCValFinalFit{cx}.sumChannelPrediction;
            end
            
            % Compute mean and concatenate predicted time series across
            % cross-val folds
            R2_crossval_mean = mean(R2_crossval,1,'omitnan');
            sumChannelPredictionCrossvalMean = mean(cat(3,...
                sumChannelPredictionCV{1},sumChannelPredictionCV{2}),3,'omitnan');
            Y_crossval_mean = mean(cat(3,...
                resultsCValInitFit{1}.Y_noOffset,resultsCValInitFit{2}.Y_noOffset),3,'omitnan');
            
            sumChannelPredictionNoOffsetCrossvalMean = mean(cat(3,...
                sumChannelPredictionCVNoOffset{1},sumChannelPredictionCVNoOffset{2}),3,'omitnan');
            Y_crossval_mean_noOffset = mean(cat(3,...
                resultsCValFinalFit{1}.Y_noOffset,resultsCValFinalFit{2}.Y_noOffset),3,'omitnan');
            
            % Chop predicted DoG time series
            [TDoG,allDoGTrials,~,~,~, catTrialsDoG] = ...
                simseq_getStimBlocksFromDataRuns(sumChannelPredictionCrossvalMean,...
                trials,er_data(1).params);
            
            % Chop simulated DoG time series
            data.predictionCatRuns = cat(1,stPred{1}.predBOLD(3:end,:),stPred{2}.predBOLD(3:end,:));
            X = data.predictionCatRuns;
            TDoG_sim = simseq_getStimBlocksFromDataRuns(X,trials,er_data(1).params);
            
            % Store results
            save(fullfile(saveFolder,sprintf('dogSimulation_%s_%s',pths.subjID, roi)),...
                'stPred','data','TDoG_sim','conditionNamesSimSeq','TDoG','allDoGTrials',...
                'B_crossval','R2_crossval_mean','Y_crossval_mean','sumChannelPredictionCrossvalMean',...
                'Y_crossval_mean_noOffset','sumChannelPredictionNoOffsetCrossvalMean','er_data','trials','params','sizeRatio')
            
        end
    end
end
  
return

