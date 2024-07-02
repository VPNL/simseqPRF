function fH = makeSupplFigure2_SuppressionDeltas(ds,lmmResults,roisToPlot, ...
                cmapROIs,subjnrs,useSTRetParams,temporalModel, spatialModel, ...
                saveFigs,saveFigDir)
% Function to reproduce supplementary figure 2: 
% panel a: Difference in suppression level (LMM  slopes) vs ROI vs stimulus condition
% panel b: Difference in suppression level (LMM  slopes) vs pRF parameters
%
% From the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University
% 
% INPUTS (required):
% - ds              : dataset
% - lmmResults      : cell (1x number of ROIs), containing a struct with fields:
%                       fixedIntercepts, fixedSlopes, 
%                       fixedIntercepts_CI, fixedSlopes_CI
% - roisToPlot      : cell with ROI names
% - cmapROIs        : color map for ROIs
% - subjnrs         : Subjects to plot
% - useSTRetParams  : (boolean) are we using supplementary spatiotemporal
%                               retinotopy data or not?
% - spatialModel    : spatial components of pRF models
% - temporalModel   : temporal components of pRF models
% - saveFigs        : save figures or not?
% - saveFigDir      : folder to save figures
%
% OUTPUTS:
% - fH         : figure handle
%
%% Check inputs
if isempty(cmapROIs) || ~exist('cmapROIs','var')
    cmapROIs = getROISummaryColors(0);
end

% Get fixed slopes and standard error
dataMat = reshape(cell2mat(lmmResults(end).fixedSlopes),4,[]);
seMat = reshape(cell2mat(lmmResults(end).fixedSlopes_SE),4,[]);

% Define inherit index based on visual hierarch
deltaLabels = {'V2-V1','V3-V2','hV4-V3','VO-hV4','V3AB-V3',...
                'IPS-V3AB','LO-V3','TO-LO'};   
deltaIdx = [2,1; ... V2 - V1
              3,2; ... V3 - V2
              4,3; ... hV4 - V3
              5,4; ... VO - hV4
              6,3; ... V3AB - V3
              7,6; ... IPS - V3AB
              8,3; ... LO - V3
              9,8]; %  TO - LO

% Compute deltas: Subtract regression slopes according to inherit index
deltasGroup = dataMat(:,deltaIdx(:,1))-dataMat(:,deltaIdx(:,2));

% Preallocate space for group 
deltasSE = {};
deltasSubj = {};

% Loop over nr of delta
for ii = 1:size(deltaIdx,1) 
    
    % Get number of contributing subjects for each area pair
    subjsArea1 = unique(ds.Subject(ismember(ds.ROI,roisToPlot(deltaIdx(ii,1)))));
    subjsArea2 = unique(ds.Subject(ismember(ds.ROI,roisToPlot(deltaIdx(ii,2)))));
    
    % preallocate space
    dataArea1 = NaN(10,4); 
    dataArea2 = dataArea1;
    
    % Get mean subject slopes and SE  
    dataArea1(subjsArea1,:) = lmmResults.subjSlopes{deltaIdx(ii,1)};
    dataArea2(subjsArea2,:) = lmmResults.subjSlopes{deltaIdx(ii,2)};
    deltasSubj{ii} = dataArea1 - dataArea2;
    deltasSE{ii} = std(deltasSubj{ii},[],1, 'omitnan')./sqrt(sum(~isnan(deltasSubj{ii}(:,1))));
end

%% Panel A: PLOT DELTA SUPPRESSION LEVELS ACROSS CONDITIONS
fH(1) = plotLMMfittedRegressionSlopesCummulative(deltasGroup,deltasSubj,deltasSE,deltaLabels,cmapROIs, saveFigs,saveFigDir);

%% Panel B: PLOT AVERAGE DELTA SUPPRESSION LEVELS ACROSS CONDITIONS AGAINST PRF SIZE/EXPONENT
% Get resampled data
output = resamplePRFParams_wReplacement(ds, roisToPlot, useSTRetParams, temporalModel, spatialModel, subjnrs);

% Plot Suppression slopes vs CST pRF size
fH(2) = plotDeltaLMMregressionSlopes_v_pRFSize(deltasGroup,deltasSubj, ...
                        output.median_resampledPRFSz, deltaIdx, cmapROIs,saveFigs,saveFigDir);
                    
%% Panel C: Plot Suppression slopes vs CST pRF exp
fH(3) = plotDeltaLMMregressionSlopes_v_pRFExp(deltasGroup,deltasSubj, ...
                        output.median_resampledCSTExp, deltaIdx, cmapROIs,saveFigs,saveFigDir);


