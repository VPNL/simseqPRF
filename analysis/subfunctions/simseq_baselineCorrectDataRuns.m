function bslCorrectedDataRuns = simseq_baselineCorrectDataRuns(dataRuns, prepostBlankTRsForBslCorrection,varargin)

[nFrames, ~] = size(dataRuns);

if nargin==3
    nFrames = varargin{1};
    nRuns = 1;
end
if nargin==4
    nFrames = varargin{1};
    nRuns = varargin{2};
end

% Take first X TRs are blank and last X TRs as blank
prepostBlank = [1:(-1*prepostBlankTRsForBslCorrection(1)), ...
    (nFrames-prepostBlankTRsForBslCorrection(2)+1):nFrames];
prepostBlank = prepostBlank + (nFrames.*[0:nRuns-1])';

% Baseline correct
bslCorrectedDataRuns = dataRuns - ...
    mean(dataRuns(prepostBlank,:),1,'omitnan');

return