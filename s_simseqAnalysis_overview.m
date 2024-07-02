%% s_simseqAnalysis_overview
% Script to rerun data analyses in the paper:
% Title:   Rethinking simultaneous suppression in visual cortex via 
%          compressive spatiotemporal population receptive fields.
% Authors: Kupers, Kim, Grill-Spector (2024).
% Journal: Nature Communications
% DOI:     XXX
%
% Requires getting MRI data from OSF (see downloadDataTableFromOSF.m)
%
% Code written by E.R. Kupers (2024) Stanford University

%% Step 0: Preprocess toonotopy data + estimate linear/CSS pRF models 
% toon_workflow.m

%% Step 1: Raw SEQ-SIM data (nifti) to mrVista gray
% simseq_preprocess_workflow_raw_to_gray

%% Step 3: Prepare data for modelfits -- make ROIS and concatenate parfiles
% s_prepareDataForModelFits

%% Step 4: Detrend, average, combine, and chop measured time series.
% s_runDataPreproc

%% Step 5: Make model predictions given pRFs
% s_runAllModelPredictions

%% Step 6: Fit model predictions to data
% s_runAllModelFits

%% Step 7: Find best R2 for exponent grid fit
% s_selectBestGridFitExp

%% Step 8: Summarize and prepare data for plotting
% s_prepareDataForPlotting.m

% For supplementary analyses:
% s_prepareSTRetDataForPlotting.m

%% Step 9: Visualize
% s_makeAllManuscriptFigures



%% Simulate difference of gaussian (DoG) pRF time series
% s_simulateDoGPRFTimeSeries.m

