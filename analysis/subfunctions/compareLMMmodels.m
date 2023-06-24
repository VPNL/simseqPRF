function compareLMMs = compareLMMmodels(modelA, modelB, modelAName,modelBName)

if ~exist('modelAName','var')
    modelAName = [];
end

if ~exist('modelBName','var')
    modelBName = [];
end

compareLMMs = compare(modelA,modelB);
fprintf('Fixed Effects comparison:\n')
fprintf('%s %s vs %s %s\n',compareLMMs.Model(1),modelAName,compareLMMs.Model(2),modelBName);
fprintf('AIC:\t %3.2f\t vs \t %3.2f\n',compareLMMs.AIC(1),compareLMMs.AIC(2))
fprintf('BIC:\t %3.2f\t vs \t %3.2f\n',compareLMMs.BIC(1),compareLMMs.BIC(2))
fprintf('Log Likelihood:\t %3.2f\t vs \t %3.2f\n',compareLMMs.LogLik(1),compareLMMs.LogLik(2))
fprintf('Likelihood ratio:\t %2.2f, \tp-value: %1.3f\n', compareLMMs.LRStat(1), compareLMMs.pValue(1))
fprintf('Likelihood ratio:\t %2.2f, \tp-value: %1.3f\n', compareLMMs.LRStat(2), compareLMMs.pValue(2))

end
