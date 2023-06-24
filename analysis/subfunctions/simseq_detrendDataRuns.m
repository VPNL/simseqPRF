function [detrendedDataRuns, params] = simseq_detrendDataRuns(dataRuns,verbose)

%% Inputs:
% DataRuns array with time series in dimension order time x runs x voxels

%% Outputs:
% detrendedDataRuns 

%% Detrend data (see also detrendTSeries as an alternative)

[nFrames,nRuns,nVox] = size(dataRuns);

% Set detrending params
params = [];
params.stim.nUniqueRep  = 1;
params.stim.nDCT        = 1; % detrending frequeny maximum (cycles per scan): 1 means 3 detrending terms, DC (0 cps), 0.5, and 1 cps
params.stim.hrfParams   = {[1.6800, 3, 2.0500], [5.4000, 5.2000, 10.8000, 7.3500, 0.3500]};
params.stim.hrfType     = 'two gammas (SPM style)';
params.stim.nFrames     = nFrames;
params.analysis.HrfMaxResponse = 1; % ??? CHECK THIS
trends = rmMakeTrends(params,verbose);

detrendedDataRuns = NaN(nFrames,nVox,nRuns);

% Do it!
for run = 1:nRuns
    currRun = squeeze(dataRuns(:,run,:));
    trendBetas = pinv(trends)*currRun;
    detrendedDataRuns(:,:,run) = currRun - trends*trendBetas;
end

return