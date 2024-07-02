function timepoints = getTimePointsTrialToAverage(selectionType, timeWindow, data)

switch selectionType
    case 'simpleOnset'
         [~,timepoints] = intersect(timeWindow,[4:12]);
        timepoints = timepoints';
    case 'variableOnset'
        [~,zeroIdx] = find(timeWindow==0);
        threshOnset = 0.1;
        blockLength = 9;
        bounds = [1:max(timeWindow)-blockLength];
        for nn = 1:size(data,1)
            dataSkipBaseline = data(nn,timeWindow>0);
            dataNorm = dataSkipBaseline+abs(dataSkipBaseline(1));
            dataNorm = dataNorm./sum(dataNorm);
            dy = cumsum(dataNorm);            
            [~,firstidx] = find(dy(bounds)>threshOnset,1,'first');
            if ~isempty(firstidx)
                onset = firstidx + zeroIdx;
            else
                onset = 4 + zeroIdx;
            end
            timepoints(nn,:) = onset:(onset - 1 + blockLength);
        end
end
