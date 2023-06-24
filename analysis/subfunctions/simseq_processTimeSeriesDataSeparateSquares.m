function [] = simseq_processTimeSeriesDataSeparateSquares(pths,hm, data, prfParams, ...
    alignMe,bslCorrect,prepostBlankTRsForBslCorrection, dt, saveDataFolder,subFolderName, saveData)
%
%   simseq_processTimeSeriesDataSeparateSquares
%
% Function to separate data of simseq experiment into separate stimulus
% square locations

%% Check voxel alignment of datasets
clear prfParamsSquare;
for h = 1:length(hm)
    if isfield(alignMe.pRFs.selected,hm{h})
        if ~exist('prfParamsSquare','var'), prfParamsSquare = prfParams; end
        fnROI = fieldnames(prfParamsSquare.(hm{h}));
        fnROI = fnROI(1:end-5);
        for ii = 1:length(fnROI)
            if strcmp(fnROI{ii},'hvolROI')
                prfParamsSquare.(hm{h}).([fnROI{ii} '_orig']) = prfParamsSquare.(hm{h}).(fnROI{ii});
            elseif strcmp(fnROI{ii},'hvolcoords')
                prfParamsSquare.(hm{h}).(fnROI{ii}) = prfParamsSquare.(hm{h}).(fnROI{ii})(:,alignMe.pRFs.selected.(hm{h}));
            else
                prfParamsSquare.(hm{h}).(fnROI{ii}) = prfParamsSquare.(hm{h}).(fnROI{ii})(alignMe.pRFs.selected.(hm{h}));
            end
        end
    end
end

% remove voxels from dataRuns and data gray coordinates
if isfield(alignMe.dataRuns, 'selected')
    if isfield(alignMe.dataRuns.selected,'lh') && ~isempty(alignMe.dataRuns.selected.lh)
        data.dataRuns.lh = data.dataRuns.lh(:,:,alignMe.dataRuns.selected.lh);
    end
    if isfield(alignMe.dataRuns.selected,'rh') && ~isempty(alignMe.dataRuns.selected.rh)
        data.dataRuns.rh = data.dataRuns.rh(:,:,alignMe.dataRuns.selected.rh);
    end
end

% Separate data by stimulus square location
for h = 1:length(hm)
    if isfield(dataRuns, hm{h})
        all = data.dataRuns.(hm{h});
        data.dataRuns.(hm{h}) = [];
        data.dataRuns.(hm{h}).all = all;
        squareidx = getIndividualSquareROIIndices(pths, prfParamsSquare, roiType,  hm{h});
        data.dataRuns.(hm{h}).square1 = data.dataRuns.(hm{h}).all(:,:,squareidx{1});
        data.dataRuns.(hm{h}).square2 = data.dataRuns.(hm{h}).all(:,:,squareidx{2});
        data.dataRuns.(hm{h}).square3 = data.dataRuns.(hm{h}).all(:,:,squareidx{3});
        data.dataRuns.(hm{h}).square4 = data.dataRuns.(hm{h}).all(:,:,squareidx{4});
        fn = fieldnames(data.dataRuns.(hm{h}));
        hemiPRFSquareIdx{h} = squareidx;
    end
end

if ~isempty(data.dataRuns.lh) && ~isempty(data.dataRuns.rh)
    data.dataRuns.both = cat(3,data.dataRuns.lh.all,data.dataRuns.rh.all);
elseif ~isempty(data.dataRuns.lh)
    data.dataRuns.both = data.dataRuns.lh.all;
elseif ~isempty(data.dataRuns.rh)
    data.dataRuns.both = data.dataRuns.rh.all;
end


if ~isempty(data.dataRuns.lh) && isempty(data.dataRuns.rh)
    %% Detrend data
    
    % Set detrending params
    params = [];
    params.stim.nUniqueRep  = 1;
    params.stim.nDCT        = 1; % detrending frequeny maximum (cycles per scan): 1 means 3 detrending terms, DC (0 cps), 0.5, and 1 cps
    params.stim.hrfParams   = {[1.6800, 3, 2.0500], [5.4000, 5.2000, 10.8000, 7.3500, 0.3500]};
    params.stim.hrfType     = 'two gammas (SPM style)';
    params.stim.nFrames     = size(data.dataRuns.both,1);
    params.analysis.HrfMaxResponse = 1; % ??? CHECK THIS
    
    [trends, ntrends, dcid] = rmMakeTrends(params,verbose);
    
    for sq = 2:length(fn)
        for run = 1:size(data.dataRuns.both,2)
            for h = 1:length(hm)
                if isfield(data.dataRuns, hm{h}) && ~isempty(data.dataRuns.(hm{h}).(fn{sq}))
                    currRun = squeeze(data.dataRuns.(hm{h}).(fn{sq})(:,run,:));
                    trendBetas = pinv(trends)*currRun;
                    data.detrendedData.(hm{h}).(fn{sq})(:,:,run) = currRun - trends*trendBetas;
                end
            end
        end
    end
    
    
    %% Concatenate corresponding unique runs
    data.dataRuns.detrendedStimRuns = [];
    for sq = 2:length(fn)
        for h = 1:length(hm)
            data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq}) = [];
            for rr = stimRun
                if isfield(data.detrendedData, hm{h}) && isfield(data.detrendedData.(hm{h}),(fn{sq}))
                    tmp = data.detrendedData.(hm{h}).(fn{sq})(:,:,pths.runOrder==rr);
                    data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq}) = cat(1,data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq}), tmp);
                end
            end
        end
    end
    
    %% Average unique runs: (1) all, (2) first half, and (3) second half
    for sq = 2:length(fn)
        for h = 1:length(hm)
            if isfield(data.dataRuns.detrendedStimRuns, hm{h}) && isfield(data.dataRuns.detrendedStimRuns.(hm{h}),(fn{sq}))
                data.meanDetrendedData.(hm{h}).(fn{sq})     = mean(data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq}),3,'omitnan');
                data.meanDetrendedDataOdd.(hm{h}).(fn{sq})  = mean(data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq})(:,:,1:2:end),3,'omitnan');
                data.meanDetrendedDataEven.(hm{h}).(fn{sq}) = mean(data.dataRuns.detrendedStimRuns.(hm{h}).(fn{sq})(:,:,2:2:end),3,'omitnan');
            end
        end
    end
    
    
    %% Baseline correct
    if bslCorrect
         nrTRsSingleRun = size(data.dataRuns.both,1);
        for sq = 2:length(fn)
            for h = 1:length(hm)
                if isfield(data.meanDetrendedData, hm{h}) && ...
                        isfield(data.meanDetrendedData.(hm{h}),(fn{sq})) && ...
                        ~isempty(data.meanDetrendedData.(hm{h}).(fn{sq}))
                    data.meanDetrendedData.(hm{h}).(fn{sq})     = simseq_baselineCorrectDataRuns(data.meanDetrendedData.(hm{h}).(fn{sq}), prepostBlankTRsForBslCorrection,nrTRsSingleRun);
                    data.meanDetrendedDataOdd.(hm{h}).(fn{sq})  = simseq_baselineCorrectDataRuns(data.meanDetrendedDataOdd.(hm{h}).(fn{sq}), prepostBlankTRsForBslCorrection,nrTRsSingleRun);
                    data.meanDetrendedDataEven.(hm{h}).(fn{sq}) = simseq_baselineCorrectDataRuns(data.meanDetrendedDataEven.(hm{h}).(fn{sq}), prepostBlankTRsForBslCorrection,nrTRsSingleRun);
                end
            end
        end
    end
    
    %% Get splithalf reliability
    for sq = 2:length(fn)
        for h = 1:length(hm)
            if isfield(data.detrendedData,hm{h}) && isfield(data.detrendedData.(hm{h}),(fn{sq}))
                if ~isempty(data.detrendedData.(hm{h}).(fn{sq}))
                    data.splitHalfRel.(hm{h}).(fn{sq}).run = simseq_getSplitHalfDataReliability( ...
                        permute(data.detrendedData.(hm{h}).(fn{sq}),[1,3,2]), 'runFlag', true, ...
                        'splitRuns', pths.runOrder);
                else
                    data.splitHalfRel.(hm{h}).(fn{sq}).run = [];
                end
            end
        end
    end
    
    %% Chop run data and predictions into trials,
    % create a new data type with average run 1,2,3 concatenated
    
    if numel(stimRun)>1
        postFix = sprintf('%d',1:max(stimRun));
        newDt = sprintf('Average_Run%s',postFix);
        cd(fullfile(simSeqSessionDir))
        d = dir(fullfile(simSeqSessionDir,'Stimuli','parfiles',sprintf('*run%s_v%d_corrected2.par',postFix,versionNr)));
        load('mrSESSION.mat');
        for dts = 1:length(dataTYPES)
            if strcmp(dataTYPES(dts).name,dt)
                motionComp_idx = dts;
            end
            if strcmp(dataTYPES(dts).name,newDt)
                dts_idx = dts;
                break;
            end
        end
        if ~exist('dts_idx','var') || isempty(dts_idx)
            dataTYPES(end+1).name = newDt;
            dts_idx = length(dataTYPES);
        end
        dataTYPES(dts_idx).scanParams(1) = dataTYPES(motionComp_idx).scanParams(1);
        dataTYPES(dts_idx).scanParams(1).parfile = fullfile('./Stimuli/parfiles',d(1).name);
        dataTYPES(dts_idx).scanParams(1).WithinScanMotion = [];
        dataTYPES(dts_idx).scanParams(1).PfileName = [];
        dataTYPES(dts_idx).blockedAnalysisParams = dataTYPES(motionComp_idx).blockedAnalysisParams(1);
        dataTYPES(dts_idx).eventAnalysisParams = dataTYPES(motionComp_idx).eventAnalysisParams(1);
        dataTYPES(dts_idx).scanParams(1).nFrames = size(data.dataRuns.both,1);
        
        save('mrSESSION.mat','mrSESSION','dataTYPES', 'vANATOMYPATH', '-append');
        
        hvol    = initHiddenGray;
        hvol    = viewSet(hvol, 'curdt', newDt);
        
        trials = er_concatParfiles(hvol, 1);
        trials.framesPerRun = size(data.dataRuns.both,1);
        
        er_data(1).params.timeWindow = repmat(totalBlockTRs, length(er_data(1).trials.condNums), 1);
        er_data(1).params.bslPeriod = repmat(baselineBlockTRs, length(er_data(1).trials.condNums), 1);
        er_data(1).params.normBsl = 0;
        er_data(1).trials = trials;
        
    else
        postFix = '1';
        er_data(1).params.timeWindow = {-2:14, -2:14, -2:14}; % BLANK, SEQ & SIM
        er_data(1).params.bslPeriod = repmat({-2:1:0}, 3, 1);
        er_data(1).params.normBsl = 0;
        trials = er_data(1).trials;
    end
    
    
    Tsquare = table();
    for h = 1:length(hm)
        for sq = 2:length(fn)
            catTrialsSquare.(hm{h}).(fn{sq}) = [];
            
            if isfield(data.meanDetrendedData, hm{h}) && isfield(data.meanDetrendedData.(hm{h}),(fn{sq}))
                if ~isempty(data.meanDetrendedData.(hm{h}).(fn{sq}))
                    
                    %% Convert struct into a table with voxel by condition, separate for data and model predictions
                    [Tall,allDataTrialsSquare,~,~,nrTrials, catTrials] = ...
                        simseq_getStimBlocksFromDataRuns(data.meanDetrendedData.(hm{h}).(fn{sq}),trials,er_data(1).params); %#ok<ASGLU>
                    
                    catTrialsSquare.(hm{h}).(fn{sq}) = catTrials;
                    Tall.Properties.VariableNames = {'Condition','DataTS','DataTSerror'};
                    
                    Tvar = table(ones(size(Tall.DataTS,1))*(h-1),ones(size(Tall.DataTS,1))*(sq-1),'VariableNames',{'Hemi_lh1_rh2','Square'});
                    
                    
                    if doCrossValidation
                        [Todd,allDataTrialsSquareOdd.(hm{h}).(fn{sq})] = ...
                            simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataOdd.(hm{h}).(fn{sq}),trials,er_data(1).params); %#ok<ASGLU>
                        Todd.Properties.VariableNames = {'Condition','DataTSOdd','DataTSOdderror'};
                        
                        [Teven,allDataTrialsSquareEven.(hm{h}).(fn{sq})] = ...
                            simseq_getStimBlocksFromDataRuns(data.meanDetrendedDataEven.(hm{h}).(fn{sq}),trials,er_data(1).params); %#ok<ASGLU>
                        Teven.Properties.VariableNames = {'Condition','DataTSEven','DataTSEvenerror'};
                        
                        % Concatenate tables
                        Tsquare = [Tvar; Tall; Todd, Teven];
                    else
                        Tsquare = [Tvar; Tall];
                    end
                    
                    if ~isempty(catTrialsSquare.(hm{h}).(fn{sq}))
                        data.splitHalfRel.(hm{h}).(fn{sq}).trial = simseq_getSplitHalfDataReliability( ...
                            catTrialsSquare.(hm{h}).(fn{sq}), 'runFlag',false,'splitRuns',[]);
                    else
                        data.splitHalfRel.(hm{h}).(fn{sq}).trial = [];
                    end
                    
                end
            end
        end
    end
   
    %% save betas
    if saveData
        
        dataToSave = {'Tsquare','allDataTrialsSquare','prfParamsSquare','hemiPRFSquareIdx','data','trials','nrTrials','catTrialsSquare'};
        if doCrossValidation
            dataToSave = {dataToSave{:}, 'allDataTrialsSquareOdd','allDataTrialsSquareEven'};
        end
        
        
        if ~exist(fullfile(saveDataFolder,subFolderName),'dir')
            mkdir(fullfile(saveDataFolder,subFolderName));
        end
        
        saveFile = fullfile(saveDataFolder,subFolderName,...
            sprintf('preprocData_%s_%s_%s_run%s_sepSquare.mat', ...
            rois{roi},roiType,spatialModel, postFix));
        
        save(saveFile, dataToSave{:},'-v7.3');
    end
    
end


return











