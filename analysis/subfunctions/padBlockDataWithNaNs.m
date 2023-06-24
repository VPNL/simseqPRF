function paddedData = padBlockDataWithNaNs(dataIn, nrTRs)

unevenTRs = cellfun(@(x) find(length(x)~=nrTRs), dataIn,'UniformOutput',false);

paddedData = dataIn;
for ii = 1:length(unevenTRs),
    if ~isempty(unevenTRs{ii}),
        diffTRs = nrTRs-length(paddedData{ii});
        paddedData{ii} = cat(1,paddedData{ii},NaN(diffTRs,1));
    end
end