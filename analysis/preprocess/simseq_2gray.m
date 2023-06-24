function simseq_2gray(sessionDir)
% simseq_2gray(sessionDir);
% 
% Function to transforms the time series from the inplane to volume anatomy 
% and averages the time series
%
% Inputs
% sessionDir    : (str) path to subject's VistaSoft session
%
% Written by EK 2021, adapted from D Finzi toonotopy code

%% Load initialized session
cd(sessionDir);
 
load mrInit_params.mat
load mrSESSION.mat

% Open hidden inplane
hi = initHiddenInplane('MotionComp_RefScan1', 1);

%% Install segmentation
% Guess class file path
classFile = fullfile(pwd, '3DAnatomy/', '*class*nii*');
w = dir(classFile);
if ~isempty(w),  
    defaultClass = fullfile(pwd, '3DAnatomy/', w(1).name);
end
if exist(defaultClass)
    installSegmentation([],[],defaultClass, 3); 
else %open GUI prompt for file name
    installSegmentation([],[],[],3); 
end

%% Transform tseries
% Now we can open a gray view
hg = initHiddenGray('MotionComp_RefScan1', 1);

% Xform time series from inplane to gray using trilinear interpolation
hg = ip2volTSeries(hi,hg,0,'linear');

%% Average tseries across duplicate runs
% averaging improves SNR by ~sqrt(runCount)

% first average all 6 functional runs
nRuns = length(mrSESSION.functionals);
% hg = averageTSeries(hg, [1:nRuns], 'Averages');

%% Concatenate all runs (to be implemented: think first about detrending runs, converting to psc, etc)

% % Create new group with single scan
% groupScans(hg, 1, 'AllRunsConcat', 'All scans concatenated')
% 
% allTSeries = [];
% for r = 1:nRuns
%     load(['./Gray/MotionComp_RefScan1/TSeries/Scan' num2str(r) '/tSeries1.mat'])
%     allTSeries = cat(3, allTSeries, tSeries);
% end
% 
% [trs, verts, ~] =  size(allTSeries);
% allTSeriesPermuted = permute(allTSeries, [3, 1, 2]);
% allTSeriesCat = reshape(allTSeriesPermuted, nRuns * trs, verts);
% 
% % Replace time series
% tSeriesPath = fullfile('./Gray/AllRunsConcat/TSeries/Scan1','tSeries.mat');
% save(tSeriesPath, 'allTSeriesCat');


%% Group identical runs
scanListRun1 = [1:2:nRuns];
scanListRun2 = [2:2:nRuns];
% scanListRun3 = [3,6];
x1 = ['Scans'];
for ll = scanListRun1;  x1 = [x1, sprintf(' %d',ll)]; end
x2 = ['Scans'];
for ll = scanListRun2;  x2 = [x2, sprintf(' %d',ll)]; end

groupScans(hg, scanListRun1, 'Run1', x1)
groupScans(hg, scanListRun2, 'Run2', x2)
% groupScans(hg, scanListRun3, 'Run3', sprintf('Scans %d %d', scanListRun3(1),scanListRun3(2)))

% replace odd runs with motion comp within/between scans
for r = scanListRun1
    load(['./Gray/MotionComp_RefScan1/TSeries/Scan' num2str(r) '/tSeries1.mat'])
    % Define new path
    tSeriesPath = fullfile('./Gray/Run1/TSeries', sprintf('Scan%d',find(r==scanListRun1)), 'tSeries.mat');
    % Replace time series
    save(tSeriesPath, 'tSeries');
end

% replace even runs with motion comp within/between scans
for r = scanListRun2
    load(['./Gray/MotionComp_RefScan1/TSeries/Scan' num2str(r) '/tSeries1.mat'])
    % Define new path
    tSeriesPath = fullfile('./Gray/Run2/TSeries', sprintf('Scan%d',find(r==scanListRun2)), 'tSeries.mat');
    % Replace time series
    save(tSeriesPath, 'tSeries');
end

% replace even runs with motion comp within/between scans
% for r = scanListRun3
%     load(['./Gray/MotionComp_RefScan1/TSeries/Scan' num2str(r) '/tSeries1.mat'])
%     % Define new path
%     tSeriesPath = fullfile('./Gray/Run3/TSeries', sprintf('Scan%d',find(r==scanListRun3)), 'tSeries.mat');
%     % Replace time series
%     save(tSeriesPath, 'tSeries');
% end

% Save session and dataTYPES
saveSession;

end

