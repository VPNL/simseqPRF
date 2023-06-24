function fitStr = simseq_LMMwrapper(dependentVar,varargin)
% Parse inputs
p = inputParser;
p.addRequired('dependentVar',@ischar);
p.addParameter('predictorVar',[], @iscell);
p.addParameter('randomInterceptVar',[], @iscell);
p.addParameter('randomSlopeVar',[], @iscell);
p.parse(dependentVar,varargin{:});

% rename variables
dependentVar         = p.Results.dependentVar;
predictorVar         = p.Results.predictorVar;
randomInterceptVar   = p.Results.randomInterceptVar;
randomSlopeVar       = p.Results.randomSlopeVar;

% Add dependent variable
fitStr = sprintf('%s ~ ',dependentVar);

% Add predictors
if exist('predictorVar','var') && ~isempty(predictorVar)
    for pp = 1:length(predictorVar)
        if pp == 1
            fitStr = [fitStr, sprintf('%s', predictorVar{pp})];
        else
            fitStr = [fitStr, sprintf(' + %s', predictorVar{pp})];
        end
    end
end

% Add random variables
if exist('randomInterceptVar','var') && ~isempty(randomInterceptVar)
    if exist('randomSlopeVar','var') && ~isempty(randomSlopeVar)
        for rr = 1:length(randomSlopeVar)
            fitStr = [fitStr, sprintf(' + (%s | %s)', randomSlopeVar{rr}, randomInterceptVar{rr})];
        end
    else
        for rr = 1:length(randomInterceptVar)
            fitStr = [fitStr, sprintf(' + (1 | %s)', randomInterceptVar{rr})];
        end
    end
end

