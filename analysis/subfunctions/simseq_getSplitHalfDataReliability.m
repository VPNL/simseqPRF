function r_splithalf = simseq_getSplitHalfDataReliability(data, varargin)
% r_splithalf = simseq_getSplitHalfReliability(data, [nIter], [runFlag])

%% Parse inputs
p = inputParser;
p.addRequired('data', @isnumeric);
p.addParameter('nIter',1000, @isnumeric);
p.addParameter('runFlag',true, @islogical);
p.addParameter('splitRuns',[],@isnumeric)
p.parse(data,varargin{:});

% rename variables
data    = p.Results.data;
nIter   = p.Results.nIter;
runFlag = p.Results.runFlag;
splitRuns = p.Results.splitRuns;

if runFlag
    % How many runs do we have?
    [nrTimePoints, nrRepeats, nrVoxels] = size(data);
    nrConditions = 1;
else
    % How many trials do we have?
    [nrTimePoints, nrRepeats, nrVoxels, nrConditions] = size(data);
end

r = NaN(nIter,nrVoxels,nrConditions);
for ii = 1:nIter
    %     for vv = 1:nrVoxels
    for cc = 1:nrConditions
        
        % Get halves by randomly selecting runs
        if isempty(splitRuns)
            order = randperm(nrRepeats);
            splitA = order(1:ceil(nrRepeats/2));
            splitB = setdiff(order,splitA);
            
            dataA = squeeze(mean(data(:,splitA,:,cc),2, 'omitnan'));
            dataB = squeeze(mean(data(:,splitB,:,cc),2, 'omitnan'));
            
        else
            dataA = [];
            dataB = [];
            for uniqueRuns = unique(splitRuns)
                order = Shuffle(find(splitRuns==uniqueRuns));
                splitA = order(1:ceil(length(order)/2));
                splitB = setdiff(order,splitA);
                
                % Average unique runs for each split
                tmpA = squeeze(mean(data(:,splitA,:,cc),2, 'omitnan'));
                tmpB = squeeze(mean(data(:,splitB,:,cc),2, 'omitnan'));
                
                % Concatenate across time
                dataA = cat(1, dataA, tmpA);
                dataB = cat(1, dataB, tmpB);
            end
        end
        % Get Pearson r correlation across time for each voxel
        r(ii,:,cc) = diag(corr(dataA,dataB));
    end
end
r_splithalf = squeeze(mean(r,1, 'omitnan'));

return


