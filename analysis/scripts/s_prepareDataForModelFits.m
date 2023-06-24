%% s_prepareDataForModelFits.m

projectDir = simseqRootPath;
subjnr     = 1;

%% Step 1: Get square ROIs
getStimROIsFromRetinotopy(projectDir,subjnr);

%% Step 2: concatenate par files with stimulus condition onset/offsets
concatParfilesUniqueRuns(projectDir, subjnr)
