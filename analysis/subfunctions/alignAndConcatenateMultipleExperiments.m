function dataRuns = alignAndConcatenateMultipleExperiments(tSeriesRun1_in,tSeriesRun2_in, coords1,coords2)
% Function to reorder voxels of fMRI time series coming from separate data
% collection sessions for the same experiment. The slice prescriptions will
% vary across session. Thus in order to average or combine data runs, we 
% want to make sure we remove non-overlapping vertices.
% This function assumes both datasets are projects onto to the same mrMesh
% and both have a set of corresponding gray coordinates.
%
% INPUTS:
% tSeriesRun1_in    : time series in (nFrames x nRuns x nVox) from
%                       experiment session 1
% tSeriesRun2_in    : time series in (nFrames x nRuns x nVox) from
%                       experiment session 2
% coords1           : mrVista Gray coords from session 1, 3 x nVoxels array
% coords2           : mrVista Gray coords from session 2, 3 x nVoxels array
%
% Written by Eline Kupers (2022) Stanford U

%% Get dimensions

[nFrames1, nRuns1, nVox1] = size(tSeriesRun1_in);
[nFrames2, nRuns2, nVox2] = size(tSeriesRun2_in);

%% Check dimensions and align coords

% Either there are more voxels in session 1..
if nVox1 > nVox2
    % Get overlap of gray coordinates with data
    [~, ai] = intersectCols(coords1,coords2);
    % Preallocate space for new tSeries array
    tSeriesRun2 = NaN(nFrames2,nRuns2,max([nVox1,nVox2]));
    % Fill in run data data
    tSeriesRun2(:,:,ai) = tSeriesRun2_in;
    % Concate runs 
    dataRuns = cat(2,tSeriesRun1_in, tSeriesRun2);

% Or in session 2..
elseif nVox2 > nVox1
    [~, ai] = intersectCols(coords2,coords1);
    tSeriesRun1 = NaN(nFrames1,nRuns1,max([nVox1,nVox2]));
    tSeriesRun1(:,:,ai) = tSeriesRun1_in;
    dataRuns = cat(2,tSeriesRun1,tSeriesRun2_in);

else % If not, we can just concatenate the runs from 2 sessions
     dataRuns = cat(2,tSeriesRun1_in,tSeriesRun2_in);
end


end