function [fLMM,fixedInterceptsFull,fixedInterceptsFull_CI,fixedInterceptsFull_SE, ...
          fixedSlopesFull, fixedSlopesFull_CI, fixedSlopesFull_SE, ...
          subjInterceptsFull, subjInterceptsFull_CI, subjInterceptsFull_SE, ...
          subjSlopesFull, subjSlopesFull_CI,subjSlopesFull_SE] = ...
            getLMMcoefficients(ds,formulaStr,fixedEffects,randomEffects,runANOVA)
        
% Check inputs for running the LMM with "effects"
% (to compare conditions against one another)       
if ~exist('runANOVA','var') || isempty(runANOVA)
    runANOVA = false;
end
        
nrSubjects   = length(unique(ds.Subject));
% nrConditions = length(unique(ds.Condition));
nrFixedVar   = length(fixedEffects);
nrRandomVar  = length(randomEffects);

%% RUN LMM
if runANOVA
    fLMM = fitlme(ds, formulaStr,'DummyVarCoding','effects');
else
    fLMM = fitlme(ds, formulaStr);
end

% Get fixed effects (fe) and random effects (re)
[fe,feNames,feStats] = fLMM.fixedEffects;
[re,reNames,reStats] = fLMM.randomEffects;
feCI = [feStats.Lower,feStats.Upper];
reCI = [reStats.Lower,reStats.Upper];
nrFixedCoefficients = length(fLMM.CoefficientNames);
nrFixedVarLevels    = size(unique(feNames),1)/(length(unique(feStats.Name)));
nrRandomGroups      = length(unique(reStats.Group));
nrRandomVarLevels   = size(unique(reNames),1)/(length(unique(reStats.Level))*nrRandomGroups)/nrRandomVar;


fixedIntercept_idx  = cellfind(feNames.Name,'(Intercept)');
hasFixedIntercept   = ~isempty(fixedIntercept_idx);

randomIntercept_idx  = cellfind(reNames.Name,'(Intercept)');
hasRandomIntercept  = ~isempty(randomIntercept_idx);

% Check for interactions
if (hasFixedIntercept & fLMM.NumPredictors > 2) || ...
        (~hasFixedIntercept & fLMM.NumPredictors > 1)
    xConditions = fLMM.NumPredictors + 1; 
    assert(xConditions==unique(fixedEffects))
    if hasFixedIntercept
        fixedIntercept_idx  = fixedIntercept_idx+[0:(xConditions-1)];
    end
    if hasRandomIntercept & (nrRandomVarLevels > 1)
        randomIntercept_idx = randomIntercept_idx+[0:(xConditions-1)];
    end
else
    xConditions = [];
end

% We only assume max one fixed and random intercept & slope & interaction factor
fixedSlope_idx      = find(~ismember([1:nrFixedCoefficients],fixedIntercept_idx));
hasFixedSlope       = ~isempty(fixedSlope_idx);

randomSlope_idx     = find(~ismember([1:size(re,1)],randomIntercept_idx(:)));
hasRandomSlope      = ~isempty(randomSlope_idx);
if hasRandomSlope
    randomSlope_idx     = reshape(randomSlope_idx',size(randomIntercept_idx,2),size(randomIntercept_idx,1))';
end
%% Get fixed effects: INTERCEPT
if hasFixedIntercept
    fixedInterceptsResidual = fe(fixedIntercept_idx)';
    fixedInterceptsResidual_CI = feCI(fixedIntercept_idx,:);
    fixedInterceptsFull_SE = feStats.SE(fixedIntercept_idx);

    if ~isempty(xConditions) && xConditions > 1
        fixedInterceptsFull    = cat(2,fixedInterceptsResidual(1),(fixedInterceptsResidual(1)+fixedInterceptsResidual(2:end)));
        fixedInterceptsFull_CI = [fixedInterceptsResidual_CI(1,:);(fixedInterceptsResidual_CI(1,:)+fixedInterceptsResidual_CI(2:end,:))];
    else
        fixedInterceptsFull    = fixedInterceptsResidual;
        fixedInterceptsFull_CI = fixedInterceptsResidual_CI;
    end
else
    fixedInterceptsFull     = [];
    fixedInterceptsFull_CI  = [];
    fixedInterceptsFull_SE  = [];
end


%% Get fixed effects: SLOPE
if hasFixedSlope
    fixedSlopesResidual    = fe(fixedSlope_idx)';
    fixedSlopesResidual_CI = feCI(fixedSlope_idx,:);
    fixedSlopesFull_SE     = feStats.SE(fixedSlope_idx);
    
    if ~isempty(xConditions)&& xConditions > 1
        fixedSlopesFull    = cat(2,fixedSlopesResidual(1),(fixedSlopesResidual(1)+fixedSlopesResidual(2:end)));
        fixedSlopesFull_CI = [fixedSlopesResidual_CI(1,:);(fixedSlopesResidual_CI(1,:)+fixedSlopesResidual_CI(2:end,:))];
    else
        fixedSlopesFull    = fixedSlopesResidual;
        fixedSlopesFull_CI = fixedSlopesResidual_CI;
    end
else
    fixedSlopesFull    = [];
    fixedSlopesFull_CI = [];
    fixedSlopesFull_SE = [];
end

%% Get random effects: 

% Random INTERCEPT 
if hasRandomIntercept
    fixedInterceptsMtx  = repmat(fixedInterceptsFull,[size(randomIntercept_idx,1),1]);
    subjInterceptsResiduals    = reStats.Estimate(randomIntercept_idx);
    subjInterceptResiduals_CI_Lower = reshape(reCI(randomIntercept_idx,1),size(randomIntercept_idx));
    subjInterceptResiduals_CI_Upper = reshape(reCI(randomIntercept_idx,2),size(randomIntercept_idx));
    subjInterceptsFull_SE = reStats.SEPred(randomIntercept_idx);
        
    % If we have single intercept per subject, as random effect
    if size(subjInterceptsResiduals,2)==1
        subjInterceptsFull = fixedInterceptsMtx+subjInterceptsResiduals; % ?
        subjInterceptsFull_CI_Lower = subjInterceptResiduals_CI_Lower; % add fixedSlopesFull_CI(1)?
        subjInterceptsFull_CI_Upper = subjInterceptResiduals_CI_Upper; % add fixedSlopesFull_CI(2)?
        subjInterceptsFull_CI = cat(3,subjInterceptsFull_CI_Lower,subjInterceptsFull_CI_Upper);

        % If we have a conditioned intercept
    elseif size(subjInterceptsResiduals,2)>1
        subjInterceptsResidual2 = cat(2,subjInterceptsResiduals(:,1), ...
            (subjInterceptsResiduals(:,2:end)+repmat(subjInterceptsResiduals(:,1),[1,(size(subjInterceptsResiduals,2)-1)])));
        subjInterceptsFull           = fixedInterceptsMtx+subjInterceptsResidual2;
        subjInterceptsFull_CI_Lower  = [subjInterceptResiduals_CI_Lower(1,:);(subjInterceptResiduals_CI_Lower(1,:)+subjInterceptResiduals_CI_Lower(2:end,:))];
        subjInterceptsFull_CI_Upper  = [subjInterceptResiduals_CI_Upper(1,:);(subjInterceptResiduals_CI_Upper(1,:)+subjInterceptResiduals_CI_Upper(2:end,:))];
        subjInterceptsFull_CI        = cat(3,subjInterceptsFull_CI_Lower,subjInterceptsFull_CI_Upper);
    else
        subjInterceptsFull_CI = cat(3,subjInterceptResiduals_CI_Lower,subjInterceptResiduals_CI_Upper);
    end
else
    subjInterceptsFull = [];
    subjInterceptsFull_CI = [];
    subjInterceptsFull_SE = [];
end

% Get random SLOPE
if hasRandomSlope
    fixedSlopesMtx      = repmat(fixedSlopesFull,[size(randomSlope_idx,1),1]);
    subjSlopesResiduals   = reStats.Estimate(randomSlope_idx);
    subjSlopesResidual_CI_Lower = reshape(reCI(randomSlope_idx,1),size(randomSlope_idx));
    subjSlopesResidual_CI_Upper = reshape(reCI(randomSlope_idx,2),size(randomSlope_idx));
    subjSlopesFull_SE    = reStats.SEPred(randomSlope_idx);

    % If we have single slope per subject, as random effect
    if size(subjSlopesResiduals,2)==1
        subjSlopesFull       = fixedSlopesMtx+subjSlopesResiduals;
        subjSlopesFull_CI_Lower = subjSlopesResidual_CI_Lower;
        subjSlopesFull_CI_Upper = subjSlopesResidual_CI_Upper;
        
    elseif size(subjSlopesResiduals,2)>1   
        subjSlopesResidual2  = cat(2,subjSlopesResiduals(:,1), ...
                                    (subjSlopesResiduals(:,2:end)+repmat(subjSlopesResiduals(:,1),1,size(randomSlope_idx,2)-1)));
        subjSlopesFull       = fixedSlopesMtx+subjSlopesResidual2;        
        subjSlopesFull_CI_Lower     = [subjSlopesResidual_CI_Lower(1,:);(subjSlopesResidual_CI_Lower(1,:)+subjSlopesResidual_CI_Lower(2:end,:))];
        subjSlopesFull_CI_Upper     = [subjSlopesResidual_CI_Upper(1,:);(subjSlopesResidual_CI_Upper(1,:)+subjSlopesResidual_CI_Upper(2:end,:))];
        subjSlopesFull_CI = cat(3,subjSlopesFull_CI_Lower,subjSlopesFull_CI_Upper);
    end
else
    subjSlopesFull = [];
    subjSlopesFull_CI = [];
    subjSlopesFull_SE = [];
end

end