function [T,allTrials,meanTrials,semTrials,nrTrials, catTrials] = ...
    simseq_getStimBlocksFromDataRuns(dataAll,trials,params)

% Chop data:              
% er_ChopTSeries2_VaryTriaLength needs data input: TRs/time by voxels
%
% VoxModel or VoxData are structs with one field for each voxel,
% each field has cells with X conditions (X-1 + blank-0),
% and there are [nr blocks/run * nr runs] +1 blocks, so all blocks
% are concatenated across unique runs (for pilot 3 it is 8*2*2)+1.
voxData = er_chopTSeries2_VaryTrialLength(dataAll,trials,params);

% Preallocate space for table
T = table; 
catTrials = [];
allTrials = {};

% Loop over conditions to get table data
for c = 1:length(voxData(1).meanTcs)-1
    condNm = c+1; % skip blank condition
    clear meanDataTrials semDataTrials 
    
    for vv = 1:length(voxData)
        theseData = voxData(vv).allTcs(:,condNm);
        theseData = theseData(~cellfun(@isempty, theseData)); % remove empty cells
        for blocknr = 1:size(theseData,1)
            nrTRs(blocknr) = length(theseData{blocknr});
        end
        if length(unique(nrTRs))==1, nrTRs=nrTRs(1);
        else,
            [~,mindx] = min(nrTRs);
            theseData{mindx} = [theseData{mindx}; NaN(diff(unique(nrTRs)),1)];
            nrTRs = size(theseData{mindx},1);
        end
        
        allTrials{c}(:,:,vv) = reshape(cell2mat(theseData), nrTRs, []);
        meanTrials(:,vv) = mean(reshape(cell2mat(theseData), nrTRs, []),2,'omitnan');
        semTrials(:,vv) = std(reshape(cell2mat(theseData), nrTRs, []),[],2,'omitnan')./sqrt(length(theseData));
        if vv ==1
            nrTrials(c) = length(theseData);
        end
    end

    % Dimensions of catTrials = time points x trials x voxels x conditions
    catTrials = cat(4,catTrials,allTrials{c});


    tmpT =  table(c,{meanTrials},{semTrials}, ...
        'VariableNames', {'Condition','TS', 'TSerror'});
    
    % Concatenate tables
    T = [T; tmpT];
end

return