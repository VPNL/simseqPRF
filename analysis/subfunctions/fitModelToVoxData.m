function results = fitModelToVoxData(X,Y,varargin)
%% Function to fit model prediction to data at single voxel level
%         params.analysis.regressionType = 'OLS';%'fracridge';
%         Y = data.meanDetrendedData;
%         X = data.predictionCatRuns;

% INPUTS:
% X                    : (int) prediction array (numTimePoints x nVox x nChannels)
% Y                    : (int) data matrix (nFrames x nVox)
% [regressionType]     : (str) Linear regression optimization type
%                          'OLS' (Ordinary Least Squares, default)
%                          'fracridge' (Fractional Ridge regression)
%                           see https://github.com/nrdg/fracridge/tree/master/matlab
% [alphaFrac]          : factor for fractional ridge regression 
% [kFolds]             : number of data folds used to crossvalidate optimal 
%                           alpha in fractional ridge regression. 
%                           THIS IS NOT FOR FINAL MODELFIT CROSSVALIDATION 
% [sumSTFlag]          : (bool) weigh Sustained and Transient channel
%                          equally and fit their sum
% [weightChannels]     : (int) 1xN vector, with N equals number of
%                          predicted BOLD channels. Using fixed weights
%                          between [0-1] to allows us to take a weighted
%                          average of all channels in an attempt to avoid
%                          poor model fits (due to collinnearity between
%                          the channels)
% [addOffsetFlag]      : (bool) add row of ones to allow for offset in GLM
% [normChannels]       : (bool) set max of predicted response to 1 or not.
%
% OUTPUTS:
% results              : struct with modelfit results, including 
%                          lm = linear model results
%                          B  = beta values
%                          R2 = Coefficient of Determination (fraction)
%                          B  = beta values
%                          sumChannelPrediction = scaled prediction by beta
%                          sumChannelPredictionNoOffset = scaled prediction
%                          by beta, but with offset removed (only if
%                          addOffsetFlag is true).
%
% Written by Eline Kupers @ Stanford U (2022)

%% Parse input variables
p = inputParser;
p.addRequired('X', @isnumeric);
p.addRequired('Y', @isnumeric);
p.addParameter('regressionType', 'OLS', @(x) any(validatestring(x,{'OLS','fracridge'})));
p.addParameter('alphaFrac',[],@isnumeric);
p.addParameter('kFolds',[],@isnumeric);
p.addParameter('sumSTFlag', false, @islogical);
p.addParameter('weightChannels', [], @isnumeric);
p.addParameter('addOffsetFlag', true, @islogical);
p.addParameter('normChannels', true, @islogical);
p.addParameter('useGPU', false, @islogical);
p.parse(X,Y,varargin{:});

% Rename variables
fnToRename = fieldnames(p.Results);
for ff = 1:length(fnToRename)
    eval([sprintf('%s = p.Results.(fnToRename{ff});', fnToRename{ff})]);
end
clear fnToRename ff

%% Prepare for modelfit

% Get dimensions of data and prediction arrays
[numTimePoints,numVoxels,numChannels] = size(X);

% Check if dimensions of data and predictions are the same
assert(isequal(size(Y,1),size(X,1)));
assert(isequal(size(Y,2),size(X,2)));

% Normalize max height across concatenated runs, for each channel
if normChannels
    X_norm = NaN(size(X));
    for chan = 1:numChannels
        X_norm(:,:,chan) = normMax(X(:,:,chan));
    end

    % Rename variable
    X = X_norm;  clear X_norm
end

% Fit sum of equal weighted Sustained and Transient channels if requested
if sumSTFlag
    X = normMax(bsxfun(@(x) plus(x,3), X));
end

% Weight channels if requested. The variable weightChannels should
% be a 1xN vector with N equal to number of channels and with
% values between 0 and 1. By forcing weights, we can sum channels
% and then fit one regressor to avoid poor fits due to collinearity.
if ~isempty(weightChannels)
    X = X*weightChannels';
    X = normMax(bsxfun(@(x) plus(x,3)));
end

% Add offset if requested
if addOffsetFlag
    X = cat(3,X,ones(numTimePoints,numVoxels,1));
end

%% Start fitting!
results = struct();

paramsToAdd = {};
if exist('alphaFrac','var') && ~isempty(alphaFrac)
    paramsToAdd([1,2]) = {'alpha',alphaFrac};
end
if ~isempty(kFolds)
    if ~isempty(paramsToAdd)
        paramsToAdd([3,4]) = {'kFolds',kFolds};
    else
        paramsToAdd([1,2]) = {'kFolds',kFolds};
    end
end

if useGPU
    X = gpuArray(X);
    Y = gpuArray(Y);
    paramsToAdd([length(paramsToAdd)+1, length(paramsToAdd)+2]) = {'useGPU',true};
end

[lm, sumChannelPrediction] = ...
    fitModelPredictionToDataWrapper(Y, X,...
    'regressionType', regressionType, ...
    paramsToAdd{:});

lmcell = struct2cell(lm);
fnBeta = find(strcmp(fieldnames(lm),'betas'));
fnR2   = find(strcmp(fieldnames(lm),'R2'));

B  = squeeze(cell2mat(lmcell(fnBeta,:,:))); %#ok<FNDSB>
R2 = squeeze(cell2mat(lmcell(fnR2,:,:))); %#ok<FNDSB>

if numVoxels > 1
    B = B';
    R2 = R2';
end

% Generate also a response without offset
if addOffsetFlag
    offsetVals = reshape(B(:,end),[1 size(B(:,end))]);
    sumChannelPredictionNoOffset = sumChannelPrediction - offsetVals.*X(:,:,end);
    Y_noOffset = Y - offsetVals;
else
    offsetVals = [];
    sumChannelPredictionNoOffset = [];
    Y_noOffset = [];
end

if useGPU
    sumChannelPrediction = gather(sumChannelPrediction);
    sumChannelPredictionNoOffset = gather(sumChannelPredictionNoOffset);
    Y_noOffset = gather(Y_noOffset);
    R2 = gather(R2);
    B = gather(B);
end

% accumulate results
results.lm = lmcell;
results.B  = B;
results.R2 = R2;
results.sumChannelPrediction = sumChannelPrediction;
results.sumChannelPredictionNoOffset = sumChannelPredictionNoOffset;
results.Y_noOffset = Y_noOffset;


end
