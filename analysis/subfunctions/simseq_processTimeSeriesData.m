function [] = simseq_processTimeSeriesData(subjnr,projectDir,varargin)
%   simseq_processTimeSeriesData(subjnr, projectDir, spatialModel, temporalModel)
%
% Function to align and chop measured fMRI data of simseq experiment
%
% INPUTS
% subjnr                            : (int) subject to analyze
% projectDir                        : (str) where does project live
% [dt]                              : (str) vistasoft dataTYPE
% [hemi]                            : (str) hemisphere to analyze, 'lh',
%                                      'rh', or 'both' (default)
% [recomputeFlag]                   : (bool) recompute (default=true) or
%                                      load cached data (false) by tc_init
%                                      voxel time series
% [preserveCoords]                  : (bool)supposed to preserve ROI
%                                      coordinates across sessions (tc init
%                                      function). I say "supposed" because
%                                      this functionality is buggy, hence
%                                      default is set to 0, and we do it
%                                      ourselves.
% [roiType]                         : (str) within ROI, one can pick a
%                                      simseq subROI. Choose from:
%                                      pilot1 'stimcorner'
%                                      pilot2 'stimcorner16'(default),'stimcorner4'
%                                      pilot3 'stimcorner4_area4sq_eccen5'
%                                             'stimcorner4_area2sq_eccen5'
%                                         'stimsquaresonly4_area4sq_eccen5'
%                                         'stimsquaresonly4_area42_eccen5'
% [verbose]                         : (bool) plot figures (default=true)
% [saveFigs]                        : (bool) save figures (default=true)
% [saveData]                        : (bool) save chopped data+fits or not
%                                       (default=true)
% [saveDataFolder]                  : (str) path to save data
% [sessionNr]                       : (str) mrVista session folder (if
%                                       there are multiple)
% [subFolder]                       : (str) subFolder for data to save
% [stimFileName]                    : (str) name of the stimulus sequence
% [doCrossValidation]               : (bool) prepare data for split-half
%                                       crossvalidation or not
%                                       (default=false)
% [extraTRsToRemoveAtBeginning]     : (int) extra trs to remove from start.
%                                       Note that we already removed TRs to
%                                       deal with scanner instability at
%                                       preprocessing step (
%                                       TRsRemovedByPreprocessing).
%                                       (default=0)
% [countDown]                       : (int) nr of TRs used for countdown,
%                                       these TRs are not part of stimulus
%                                       sequence, but are part of data run
%                                       (default=6 TRs).
% [TRsRemovedByPreprocessing]       : (int) nr of TRs removed for each data
%                                       run as part of preprocessing
%                                       (default=8 TRs)
% [prepostBlankTRsForBslCorrection] : (int) TRs at the beginning and end of
%                                       a run used to correct for baseline
% [totalBlockTRs]                   : (cell) TRs before/after onset of
%                                       stimulus block (trials.onsetFrames)
%                                       used to chop runs (default={-4:1:18})
% [baselineBlockTRs]                : (cell) TRs before onset of stimulus
%                                        block used to define baseline
%                                        before block onset (default={-4:1:0})
% [separateSquareROIData]           : (bool) if we want to average run data
%                                       and trial data across the pRFs that
%                                       fall within one stimulus square.
%                                       (default=false)
% [spatialpRFModel]                 : (str) what spatial pRF model results
%                                       to load? Choose: 'cssFit' (default)
%                                       or 'onegaussianFit'.
% [veThreshpRF]                     : (int) lower variance explained
%                                       threshold for pRFs , as fraction
%                                       between [0,1] (default=0.1)
% [stimHz]                          : (int) sample rate of stimulus sequence
%                                       in Hz (default=60Hz).

%% 1. Set params and folders
p = inputParser;
p.addRequired('subjnr', @isnumeric);
p.addRequired('projectDir',@ischar);
p.addParameter('dt','MotionComp_RefScan1', @ischar);
p.addParameter('hemi','both', @(x) any(validatestring(x,{'lh','rh','both'})));
p.addParameter('recomputeFlag', true, @islogical);
p.addParameter('preserveCoords', false, @islogical);
p.addParameter('roiType','stimcorner16',@ischar);
p.addParameter('roiIdx','all', @(x) (ischar(x) || isnumeric(x) || isscalar(x) || isvector(x)));
p.addParameter('bslCorrect',false,@islogical);
p.addParameter('verbose', true, @islogical);
p.addParameter('saveData', true, @islogical);
p.addParameter('saveDataFolder',[], @ischar);
p.addParameter('sessionNr',2, @isnumeric);
p.addParameter('subFolder',[], @ischar);
p.addParameter('stimFileName',[], @ischar);
p.addParameter('doCrossValidation',false,@islogical);
p.addParameter('extraTRsToRemoveAtBeginning',0,@isnumeric);
p.addParameter('countDown',6,@isnumeric);
p.addParameter('TRsRemovedByPreprocessing',8,@isnumeric);
p.addParameter('prepostBlankTRsForBslCorrection',[-10,5],@isnumeric);
p.addParameter('totalBlockTRs',{-4:1:18},@iscell);
p.addParameter('baselineBlockTRs',{-4:1:0},@iscell);
p.addParameter('separateSquareROIData',false,@islogical);
p.addParameter('spatialpRFModel','cssFit',@ischar);
p.addParameter('veThreshpRF', 0.1, @isnumeric);
p.addParameter('stimHz',60,@isnumeric);
p.parse(subjnr,projectDir,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff 

%% Set subject params
cd(fullfile(projectDir,'data','simseq/'))
pths    = getSubjectPaths(projectDir,subjnr);
stimRun = unique(pths.runOrder);
versionNr  = pths.expversionNr;  % subj09 has v1, subj04 has v2, subj01 has v3

% Define folders to save beta values, data and figures
if saveData && isempty(saveDataFolder)
    saveDataFolder = fullfile(pths.simseqResultsDir,'preprocData');
end

% Get hemis
if strcmp(hemi, 'both'); hm = {'lh','rh'}; else hm = {hemi}; end

% Get ROIs depending on selected hemi(s)
if strcmp(hemi,'both')
    rois       = pths.definedROIsBOTH;
elseif strcmp(hemi,'lh')
    rois       = pths.definedROIsLH;
elseif strcmp(hemi,'rh')
    rois       = pths.definedROIsRH;
end

% Check for selected ROIs
if ~strcmp(roiIdx,'all') && isnumeric(roiIdx)
    rois = rois(roiIdx);
elseif ~strcmp(roiIdx,'all') && ischar(roiIdx)
    rois = find(ismember(rois,{roiIdx}));
end

%% Loop over selected rois
for roi = 1:length(rois) % find(ismember(rois,'VO1')) % %
    % Define subfolder if not defined already
    if ~isempty(subFolder)
        subFolderName = sprintf('%s%d_v%d',subFolder, max(stimRun), pths.expversionNr);
    end
    % clear ROI data struct
    clear data
    data.dataRuns.lh = [];
    data.dataRuns.rh = [];
    
    % Get simseq data from all runs, per hemi
    for h = 1:length(hm)
        thisHemi = (cell2mat(hm(h)));
        
        % Get data from both sessions
        for ses = pths.sessionNrs
            if ses > 1 % update paths if there are more than 1 scan session
                pths       = getSubjectPaths(projectDir,subjnr,ses);
            end
            % Define session dir
            simSeqSessionDir = fullfile(pths.dataDirSimSeq, pths.subjID, pths.session);
            cd(simSeqSessionDir);
            
            % Define roi name
            roiName = sprintf('%s_simseq%s',pths.preferredRoiName(thisHemi,rois{roi}),roiType);
            
            % Get event related data struct
            tmp_data = simseq_loadEventRelatedData(simSeqSessionDir, roiName, ...
                'dataType',dt, 'scansToLoad',pths.scansToLoad{ses}, ...
                'recomputeFlag', recomputeFlag, 'preserveCoordsFlag', preserveCoords);
            
            if ~isempty(tmp_data)
                er_data(ses==pths.sessionNrs) = tmp_data;
                
                % Remove vox data for now, we will recompute this later with
                % appropriate trial lengths
                er_data(ses==pths.sessionNrs).voxData = [];
                
                % Update total scans from current scan session
                totalScans  = length(pths.scansToLoad{ses});
                
                % Reshape tSeries into separate runs
                tSeriesRun{ses==pths.sessionNrs} = reshape(er_data(ses==pths.sessionNrs).tSeries, ...
                    size(er_data(ses==pths.sessionNrs).tSeries,1)/totalScans, ...
                    totalScans, []);
                voxCoords{ses==pths.sessionNrs,h} = er_data(ses==pths.sessionNrs).coords;
            else
                tSeriesRun{ses==pths.sessionNrs} = [];
                voxCoords{ses==pths.sessionNrs,h} = [];
            end
        end
        
        % Check differences in nr voxels if there are two scan sessions
        if length(pths.sessionNrs)>1
            data.dataRuns.(thisHemi) = ...
                alignAndConcatenateMultipleExperiments(tSeriesRun{1},tSeriesRun{2}, er_data(1).coords,er_data(2).coords);
            
        else % If no differences, then just add the time series of this run
            data.dataRuns.(thisHemi) = tSeriesRun{ses==pths.sessionNrs};
        end
    end
    
    % If there is no data within the ROI, we skip
    if ~isempty(tSeriesRun{:})
        
        % Free some memory
        clear tSeriesRun;
        
        % Store data gray coordinates
        coords.dataRuns.lh = voxCoords{1,1};
        if size(voxCoords,2)>1, coords.dataRuns.rh = voxCoords{1,2}; end
        
        %% Get pRFs
        roiNameNoHemi = sprintf('%s_toon_simseq%s',rois{roi},roiType);
        [prfParams, ~] = simseq_loadPRFParams(pths, 'Averages', hemi, ...
            roiNameNoHemi, spatialpRFModel, veThreshpRF, [], 'simseqROIs');
        
        %% Get pRF gray coordinates
        if isfield(prfParams, 'lh')
            coords.pRFs.lh = prfParams.lh.hvolROI.coords;
            selectedpRFs.predVEmask.lh = prfParams.lh.veMask;
        end
        if isfield(prfParams, 'rh')
            coords.pRFs.rh = prfParams.rh.hvolROI.coords;
            selectedpRFs.predVEmask.rh = prfParams.rh.veMask;
        end
        
        %% If requested remove even more TRs at the start due to scanner instability
        % we already removed 8 during preprocessing
        if extraTRsToRemoveAtBeginning>0
            if ~isempty(data.dataRuns.lh)
                data.dataRuns.lh(1:extraTRsToRemoveAtBeginning,:,:) = [];
            end
            if ~isempty(data.dataRuns.rh)
                data.dataRuns.rh(1:extraTRsToRemoveAtBeginning,:,:) = [];
            end
        end
        
        %% Check if data in scan has more TRs than the stimulus, if so, we trim
        % TRs at the end. If we don't do this, the modelpredictions won't align
        % with the data when concatenating unique runs.
        cd(simSeqSessionDir)
        if ~isempty(stimFileName)
            if numel(stimRun)>1
                stim = load(sprintf('./Stimuli/%s1_v%d.mat', stimFileName, pths.expversionNr));
                % Gather stimulus sequence per condition, usefull for plotting later
                for cond = 1:length(stim.params.sequenceAllTruncated)
                    stimMS{cond} = squeeze(stim.params.sequenceAllTruncated{cond}(1,1,:));
                end
            else
                stim = load(sprintf('./Stimuli/%s_run%d.mat',stimFileName,stimRun));
                stimMS{1} = stim.sequence;
            end
        else, stimMS = {};
        end
        
        %% Use stimulus sequence to check for extra TRs at the end
        stimLength = (size(stim.images,3)/stimHz) + countDown - TRsRemovedByPreprocessing;
        if ~isempty(data.dataRuns.lh)
            extraTRsToRemoveAtEnd = size(data.dataRuns.lh,1) - stimLength;
            toDelete = size(data.dataRuns.lh,1):-1:(size(data.dataRuns.lh,1)-extraTRsToRemoveAtEnd+1);
            data.dataRuns.lh(toDelete,:,:)=[];
            assert(stimLength==size(data.dataRuns.lh,1))
        end
        if ~isempty(data.dataRuns.rh)
            extraTRsToRemoveAtEnd = size(data.dataRuns.rh,1) - stimLength;
            toDelete = size(data.dataRuns.rh,1):-1:(size(data.dataRuns.rh,1)-extraTRsToRemoveAtEnd+1);
            data.dataRuns.rh(toDelete,:,:)=[];
            assert(stimLength==size(data.dataRuns.rh,1))
        end
        % Free some memory
        clear stim stimLength toDelete
        
        %% Check voxel alignment of datasets
        alignMe = simseq_checkVoxelCorrespondance(coords, selectedpRFs);
        
        % remove voxels from dataRuns and data gray coordinates
        if isfield(alignMe.dataRuns, 'selected')
            if isfield(alignMe.dataRuns.selected,'lh') && ~isempty(alignMe.dataRuns.selected.lh)
                data.dataRuns.lh = data.dataRuns.lh(:,:,alignMe.dataRuns.selected.lh);
            end
            if isfield(alignMe.dataRuns.selected,'rh') && ~isempty(alignMe.dataRuns.selected.rh)
                data.dataRuns.rh = data.dataRuns.rh(:,:,alignMe.dataRuns.selected.rh);
            end
        end
        
        %% If we want separate data by square ROIs, go here:
        if separateSquareROIData
            simseq_processTimeSeriesDataSeparateSquares(pths, hm, data, ...
                prfParams, alignMe,bslCorrect,prepostBlankTRsForBslCorrection, ...
                dt, saveDataFolder,subFolderName, saveData);
        end

        %% Concatenate left and right hemi into both hemis
        % Run data = time by runs by voxels (in 1s TRs)
        if ~isempty(data.dataRuns.lh) && ~isempty(data.dataRuns.rh)
            data.dataRuns.both = cat(3,data.dataRuns.lh,data.dataRuns.rh);
        elseif ~isempty(data.dataRuns.lh)
            data.dataRuns.both = data.dataRuns.lh;
        elseif ~isempty(data.dataRuns.rh)
            data.dataRuns.both = data.dataRuns.rh;
        end
        
        % Check if dataRuns hemi fields are empty or not
        if isempty(data.dataRuns.lh) && isempty(data.dataRuns.rh)
            break;
        else
            
            %% Detrend data (see also detrendTSeries as an alternative)
            [detrendedDataRuns,detrendParams] = simseq_detrendDataRuns(data.dataRuns.both,verbose);
            data.detrendedData.both = detrendedDataRuns;
            clear detrendedDataRuns
            
            %% Concatenate corresponding unique runs
            data.dataRuns.detrendedStimRuns = [];
            if numel(pths.runOrder)==1  && size(data.dataRuns.both,2)>1
                pths.runOrder = ones(1,size(data.dataRuns.both,2));
            end
            
            for rr = stimRun
                tmp = data.detrendedData.both(:,:,pths.runOrder==rr);
                data.dataRuns.detrendedStimRuns = cat(1,data.dataRuns.detrendedStimRuns, tmp);
            end
            
            %% Average unique runs: (1) all, (2) first half, and (3) second half
            data.meanDetrendedData     = mean(data.dataRuns.detrendedStimRuns,3,'omitnan');
            data.meanDetrendedDataOdd  = mean(data.dataRuns.detrendedStimRuns(:,:,1:2:end),3,'omitnan');
            data.meanDetrendedDataEven = mean(data.dataRuns.detrendedStimRuns(:,:,2:2:end),3,'omitnan');
            
            %% Baseline correct (for now set it to false, because we remove baseline at fitting stage
            if bslCorrect
                nrTRsSingleRun = size(data.dataRuns.both,1);
                data.meanDetrendedData = simseq_baselineCorrectDataRuns(data.meanDetrendedData, ...
                    prepostBlankTRsForBslCorrection, nrTRsSingleRun);
                
                data.meanDetrendedDataOdd = simseq_baselineCorrectDataRuns(data.meanDetrendedDataOdd, ...
                    prepostBlankTRsForBslCorrection,nrTRsSingleRun);
                
                data.meanDetrendedDataEven = simseq_baselineCorrectDataRuns(data.meanDetrendedDataEven, ...
                    prepostBlankTRsForBslCorrection,nrTRsSingleRun);

            end
            
            %% Get split half reliability across runs
            data.splitHalfRel.both.run = simseq_getSplitHalfDataReliability( ...
                permute(data.detrendedData.both,[1,3,2]), 'runFlag', true, ...
                'splitRuns', pths.runOrder);
            
            data.splitHalfRel.both.runCV = simseq_getSplitHalfDataReliability( ...
                permute(data.dataRuns.detrendedStimRuns,[1,3,2]), ...
                'runFlag',true, 'splitRuns', []);
            
            %% Chop run data and predictions into trials,
            % create a new data type with average run 1,2,3 concatenated
            dataAll = data.meanDetrendedData;
            dataOdd = data.meanDetrendedDataOdd;
            dataEven = data.meanDetrendedDataEven;
            
            if numel(stimRun)>1
                postFix = sprintf('%d',1:max(stimRun));
                newDt = sprintf('Average_Run%s',postFix);
                cd(fullfile(simSeqSessionDir))
                d = dir(fullfile(simSeqSessionDir,'Stimuli','parfiles',sprintf('*run%s_v%d_corrected2.par',postFix,versionNr)));
                load('mrSESSION.mat');
                %             global dataTYPES;
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
                dataTYPES(dts_idx).scanParams(1).nFrames = size(dataAll,1);
                
                save('mrSESSION.mat','mrSESSION','dataTYPES', 'vANATOMYPATH', '-append');
                
                hvol    = initHiddenGray;
                hvol    = viewSet(hvol, 'curdt', newDt);
                
                trials = er_concatParfiles(hvol, 1);
                trials.framesPerRun = size(dataAll,1);
                
                er_data(1).params.timeWindow = repmat(totalBlockTRs, length(er_data(1).trials.condNums), 1);
                er_data(1).params.bslPeriod = repmat(baselineBlockTRs, length(er_data(1).trials.condNums), 1);
                er_data(1).params.normBsl = 0;
                er_data(1).trials = trials;
            else
                postFix = '1';
                er_data(1).params.timeWindow = {-2:14, -2:14, -2:14}; % BLANK, SEQ & SIM
                er_data(1).params.bslPeriod = repmat({-2:1:0}, 3, 1);
                er_data(1).params.normBsl = 0;
                er_data(1).trials.onsetFrames    = er_data(1).trials.onsetFrames - removedTRsAtStart;
                er_data(1).trials.onsetFrames(1) = 1;
                er_data(1).trials.onsetSecs      = er_data(1).trials.onsetFrames - removedTRsAtStart;
                er_data(1).trials.onsetSecs(1)   = 0;
                er_data(1).trials.framesPerRun   = er_data(1).trials.framesPerRun - removedTRsAtStart;
                trials = er_data(1).trials;
            end
            
            %% Convert struct into a table with voxel by condition, separate for data and model predictions
            [Tall,allDataTrials,~,~,nrTrials, catTrials] = ...
                simseq_getStimBlocksFromDataRuns(dataAll,trials,er_data(1).params); %#ok<ASGLU>
            Tall.Properties.VariableNames = {'Condition','DataTS','DataTSerror'};
            
            if doCrossValidation
                [Todd,allDataTrialsOdd] = ...
                    simseq_getStimBlocksFromDataRuns(dataOdd,trials,er_data(1).params); %#ok<ASGLU>
                Todd.Properties.VariableNames = {'Condition','DataTSOdd','DataTSOdderror'};
                [Teven,allDataTrialsEven] = ...
                    simseq_getStimBlocksFromDataRuns(dataEven,trials,er_data(1).params); %#ok<ASGLU>
                Teven.Properties.VariableNames = {'Condition','DataTSEven','DataTSEvenerror'};
                
                % Concatenate tables
                T = [Tall; Todd, Teven];
            else
                T = Tall;
            end

            % Compute splithalf reliability across trials
            data.splitHalfRel.both.trial = simseq_getSplitHalfDataReliability( ...
                catTrials, 'runFlag',false,'splitRuns',[]);
            
            %% save betas
            if saveData
                dataToSave = {'T','allDataTrials','data',...
                    'nrTrials','er_data','stimMS','alignMe',...
                    'prfParams','detrendParams','trials'};
                
                if doCrossValidation
                    dataToSave = {dataToSave{:},'allDataTrialsOdd','allDataTrialsEven'};
                end
                
                if ~exist(fullfile(saveDataFolder,subFolderName),'dir')
                    mkdir(fullfile(saveDataFolder,subFolderName));
                end
                
                saveFile = fullfile(saveDataFolder,subFolderName,...
                    sprintf('preprocData_%s_%s_%s_run%s.mat', ...
                    rois{roi},roiType,spatialpRFModel, postFix));
                
                save(saveFile, dataToSave{:},'-v7.3');
            end
            
        end
    end
end

return

