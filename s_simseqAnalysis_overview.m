%% s_simseqAnalysis_overview
%
% Requires getting MRI data from flywheel, 
% renaming files using "renameFiles_flywheel2oak"
% Preprocess toonotopy data + pRF model using toon_workflow.m

%% Step 1: Raw data (nifti) to mrVista gray
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

%% Step 9: Visualize
% s_makeAllManuscriptFigures