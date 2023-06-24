function anal = er_chopTSeries2_VaryTrialLength(tSeries,trials,params,varargin);
% anal = er_chopTSeries2(tSeries,trials,params,[options]);
%
% View-independent verson of er_chopTSeries.
% chops up entered tSeries according to the assigned parfiles.
%
% Inputs are tSeries (matrix where rows are time points, cols are
% different tSeries, say from diff't rois/voxels), and trials
% (a struct obtained from er_concatParfiles, containing design matrix
% information).
%
% returns an analysis struct with the following fields:
%
%   wholeTc: vector of whole time course, concatenated across scans,
%            for selected scans, voxels
%   allTcs:  3D matrix of time courses for every trial. The rows
%            are different time points, the columns are different
%            trials, and the slices (z dim) are different conditions.
%            As with other trials, the data is taken from the specified
%            time window, incl. prestim frames.
%   meanTcs: matrix of mean time courses for each condition. Rows 
%            are different time points, columns are different conds.
%   sems:    corresponding standard errors of the mean for meanTcs.
%   timeWindow: vector specifying the time in seconds, relative to the
%               trial onset, from which each trial / mean trial time 
%               course is taken. [default is -4:16].
%   peakPeriod: time, in seconds, where the HRF is expected to peak.
%               This will be used for t-test and amplitude results below.
%               [default is 8:14].
%   bslPeriod:  time, in seconds, to consider as a baseline period. 
%               [default is 2:6].
%   amps:       Estimated amplitudes for each trial, taken
%               as the mean amplitude during the peak period
%               minus the mean amplitude during the baseline period.
%               This is a 2D matrix where rows are different trials
%               and columns are different conditions.
%   relamps:    relative fMRI amplitudes (e.g. dot products with the
%               mean time course) for each trial, in the same format
%               as amps. These should be less sensitive to the choice
%               of peak and baseline periods than amps, but may not be
%               accurate if different conditions have fundamentally
%               different response shapes (e.g., if there's a baseline
%               period where everything is decreasing and all other 
%               conditions are increasing.)
%   Hs:         1 x nConds binary vector reflecting whether
%               each condition had a mean significant activation
%               during peak relative to baseline periods
%               (one-sided t-test alpha = 0.05 by default but
%               may be entered as an optional argument).
%   ps:         Corresponding p-values for the t-tests for Hs.
%   SNR:        Mean signal to noise ratio for all trials.
%               (mean signal for peakPeriod - bslPeriod)/(stdev bslPeriod)
%   SNRdb:      Expression of SNR as decibels: 20 * log10(SNR).
%
% Many params can be entered as optional arguments. In these cases,
% call them with the form ...'arg name',[value],.... The fields
% that can be entered are: timeWindow, peakPeriod, bslPeriod, alpha.
%
% Further options:
%   barebones:          do only a minimal analysis, extracting
%                       mean time courses, amplitudes, and SEMs.
%                       This is useful for across-voxel analyses.
%
%   normBsl,[1 or 0]:   if 1, will align trials during the baseline
%                       period [default is 1].
%   alpha,[val]:        alpha value for t-tests [default 0.05].
%   onsetDelta,[val]:   automatically shift the onsets in parfiles
%                       relative to the time course by this amount
%                       (e.g. to compensate for HRF rise time).
%   'mrvWaitbar':          put up a mrvWaitbar instead of showing progress in
%                       the command line.
%   'findPeaks':        when calculating response amplitudes, figure
%                       out the peak amplitude separately for each 
%                       condition, taking the peak and surrounding 2 points.
%   [params struct]:    input a struct containing all of these
%                       params. See er_getParams.
%
% 06/17/04 ras: wrote it.
% 07/28/04 ras: clarified annotation.
% 01/25/05 ras: can now input params struct.
% 05/25/05 ras: fixed calculation of relative amplitudes.
if ieNotDefined('params')
    params = er_defaultParams;
end

%%%%% params/defaults %%%%%
barebones = 1;                  % if 0, do full analysis; if 1, do minimal analysis
normBsl = params.normBsl;       % flag to zero baseline or not
alpha = params.alpha;           % threshold for significant activations
bslPeriod = params.bslPeriod;   % period to use as baseline in t-tests, in seconds
peakPeriod = params.peakPeriod; % period to look for peaks in t-tests, in seconds
timeWindow = params.timeWindow; % seconds relative to trial onset to take for each trial
onsetDelta = params.onsetDelta; % # secs to shift onsets in parfiles, relative to time course
snrConds = params.snrConds;     % For calculating SNR, which conditions to use (if empty, use all)
waitbarFlag = 0;        % flag to show a graphical mrvWaitbar to show load progress
findPeaksFlag = 0;      % when computing amps, find peak period separately for each cond
TR = trials.TR;

%%%%% parse the options %%%%%
varargin = unNestCell(varargin);
for ii = 1:length(varargin)
    if isstruct(varargin{ii})
        % assume it's a params struct
        names = fieldnames(varargin{ii});
        for jj = 1:length(names)
            cmd = sprintf('%s = varargin{i}.%s;',names{jj},names{jj});
            eval(cmd);
        end
            
    elseif ischar(varargin{ii})
        switch lower(varargin{ii})
        case 'barebones', barebones = 1;
        case 'normbsl', normBsl = varargin{ii+1};
        case 'alpha', alpha = varargin{ii+1};
        case 'peakperiod', peakPeriod = varargin{ii+1};
        case 'bslperiod', bslPeriod = varargin{ii+1};
        case 'timewindow', timeWindow = varargin{ii+1};
        case 'onsetdelta', onsetDelta = varargin{ii+1};
        case 'snrconds', snrConds = varargin{ii+1};
        case 'waitbarflag', waitbarFlag = 1;
        case 'findpeaks', findPeaksFlag = 1;
        case 'addtrstotimewindow', addTRsToTimeWindow = varargin{ii+1};
        otherwise, % ignore
        end
    end
end

%%%%% account for format of time series
% in this code compared to elsewhere (dumb, I know, but there's a reason)
if size(tSeries,2)==1 & size(tSeries,1) > 1 % column vector
    wholeTc = tSeries'; % flip the col vector around
elseif (size(tSeries,1)>1) & (size(tSeries,2)>1)
    % what the heck, let's go recursively on the columns
    for tt = 1:size(tSeries,2)
        anal(tt) = er_chopTSeries2_VaryTrialLength(tSeries(:,tt),trials,params,varargin);
    end
    return
else
    wholeTc = tSeries;
end
clear tSeries;

% account for onset shift
if mod(onsetDelta,TR) ~= 0
    % ensure we shift by an integer # of frames
    onsetDelta = TR * round(onsetDelta/TR);
end
trials.onsetSecs = trials.onsetSecs + onsetDelta;
trials.onsetFrames = trials.onsetFrames + onsetDelta/TR;

%%%%% get nConds from trials struct
condNums = unique(trials.cond(trials.cond >= 0)); 
nConds = length(condNums);

%%%%% get a set of label names, if they were specified in the parfiles
for ii = 1:nConds
    ind = find(trials.cond==condNums(ii));
    labels{ii} = trials.label{ind(1)};
end

for ii = 1:length(trials.cond)-1
        cond(ii) = trials.cond(ii);
end


t1 = []; t2 = [];
for ii = 1:nConds
    t1(ii) = min(params.bslPeriod{ii}); 
    bslPeriod_min(ii) = t1(ii);
    bslPeriod_max(ii) = max(params.bslPeriod{ii}); 

    t2(ii) = max(params.timeWindow{ii}); 
end

%%%%% convert params expressed in secs into frames
% first, find the min and max time in seconds expressed in time window:
% with this change, we ignore the intervening time-window values (e.g.,
% timeWindow=[-4 20] and timeWindow=[-4:2:20] produce the same result, so
% we're more flexible when it comes to non-integer TRs)
% t1 = min(timeWindow);  t2 = max(timeWindow);

% second, re-express t1 and t2 in terms of TRs / MR frames
f1 = fix(t1 / TR);  f2 = fix(t2 / TR);


% figure out auxilliary values: how many frames before the stimulus starts?
prestim = max( -f1, 0 ); % if f1 > 0, there is no prestim (only zuul :)

% which frames reflect the peak and baseline periods?
pk1 = fix( min(peakPeriod) / TR );  
pk2 = fix( max(peakPeriod) / TR );

bs1 = fix( bslPeriod_min / TR );
bs2 = fix( bslPeriod_max / TR ); 

% now have a 'frame window' specifying the time window in integer MR frames
for ii = 1:nConds
    frameWindow{ii} = f1(ii):f2(ii);
    peakFrames{ii} = find(ismember(frameWindow{ii}, pk1:pk2)); 
    bslFrames{ii} = find(ismember(frameWindow{ii}, bs1(ii):bs2(ii))); 
end

%%%%% remove trials at the very end of a scan, without
%%%%% enough data to fill the time window
cutOff = find(trials.onsetFrames+frameWindow{1}(end) > length(wholeTc)+1);
if ~isempty(cutOff)
    keep = setdiff(1:length(trials.cond),cutOff);
    trials.cond = trials.cond(keep);
    trials.onsetFrames = trials.onsetFrames(keep);
    trials.onsetSecs = trials.onsetSecs(keep);
end

%%%%% build allTcs matrix of time points x trials x  conditions
%%%%% take (frameWindow) secs from each trial
allTcs = {};

for ii = 1:nConds
   cond = condNums(ii);
   ind = find(trials.cond==cond);
   for jj = 1:length(ind)
       tstart = trials.onsetFrames(ind(jj));
       tend = min([tstart+frameWindow{ii}(end),length(wholeTc)]);
       trialTimePoints = tstart:tend;

       % add prestim
       if tstart < prestim(ii)+1
           % for 1st trial, no baseline available -- set to 0
           allTcs{jj,ii} = [zeros(1,prestim(ii)) wholeTc(trialTimePoints)]';
       else
           % augment the range by previous [prestim] frames
           fullrng = trialTimePoints(1)-prestim:trialTimePoints(end);
           allTcs{jj,ii} = wholeTc(fullrng)';   
       end
       
       % remove baseline estimate, if selected
       if normBsl
           % estimate DC offset by prestim baseline vals
           DC = nanmean(allTcs{jj,ii}(bslFrames{ii},:));
           allTcs{jj,ii} = allTcs{jj,ii} - DC;
       end
   end 
end 

%%%%% find 'empty' trials, set to NaNs
% (Empty trials will result if some conditions have more
% trials than others -- in the conditions w/ fewer trials,
% the allTcs matrix will be padded with 0s to keep it a cube).
for y = 1:size(allTcs,1)
    for z = 1:size(allTcs,2)
        if all(allTcs{y,z}==0)
            allTcs{y,z} = NaN(size(allTcs{y,z}));
        end
    end
end

%%%%% get mean time courses, sems for each condition
meanTcs = cell(1,nConds);
sems = cell(1,nConds);
maxNTrials = size(allTcs,1);
trialData = NaN(size(allTcs,1),size(allTcs{1,1},1));

for ii = 1:nConds
    clear trialData;
    for mm=1:size(allTcs,1)
        emptryTrial(mm) = isempty(allTcs{mm,ii});
        if ~emptryTrial(mm)
            trialData(mm,1:size(allTcs{mm,ii},1)) = allTcs{mm,ii}';
        end
    end
    nTrials(ii) = sum(~emptryTrial);
    if nTrials(ii) > 1
        meanTcs{ii} = nanmean(trialData)';
        sems{ii} = nanstd(trialData)' ./ sqrt(nTrials);
    else
        meanTcs{ii} = trialData';
    end
end

% %%%%% calc amplitudes, do t-tests of post-baseline v. baseline
for ii = 1:nConds
    if findPeaksFlag==1
        % find the peak separately for
        % each condition
        maxVal = max(meanTcs{ii}(2:end));
        maxT = find(meanTcs{ii}==maxVal);
        for jj = 1:nTrials(ii)
            peak{jj,ii} = allTcs{jj,ii}(maxT-1:maxT+1);
            bsl{jj,ii}  = allTcs{jj,ii}(bslFrames);
            if isempty( peak{jj,ii} ) || isempty( bsl{jj,ii} )
                amps{jj,ii} = NaN;
            else
                amps{jj,ii} = (mean(peak{jj,ii}) - mean(bsl{jj,ii}))';
            end
        end
    else
        for jj = 1:nTrials(ii)
            peak{jj,ii} = allTcs{jj,ii}(peakFrames{ii});
            bsl{jj,ii} = allTcs{jj,ii}(bslFrames{ii});
            if isempty( peak{jj,ii} ) || isempty( bsl{jj,ii} )
                amps{jj,ii} = NaN;
            else
                amps{jj,ii} = (mean(peak{jj,ii}) - mean(bsl{jj,ii}))';
            end
        end
    end

end

%%%%% assign everything to the output struct
anal.wholeTc = double(wholeTc);
anal.allTcs = allTcs;
anal.meanTcs = meanTcs;
anal.sems = sems;
anal.labels = labels;
anal.timeWindow = timeWindow;
anal.peakPeriod = peakPeriod;
anal.bslPeriod = bslPeriod;
anal.condNums = condNums;
anal.amps = amps;


return
