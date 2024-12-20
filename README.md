# simseqPRF
Code accompanying the paper "Rethinking simultaneous suppression in visual cortex via compressive spatiotemporal population receptive fields" by Kupers, Kim, & Grill-Spector (2024) published in Nature Communications.

## Goal
The goal of this project is to operationalize and elucidate how simultaneous and sequential visual stimuli are processed within population receptive fields, in space and time, and generate a lower response for simultaneous over sequential stimulus presentations of otherwise similar colorful peripheral square stimuli.

## Paper citation
Kupers, E.R., Kim, I. & Grill-Spector, K. Rethinking simultaneous suppression in visual cortex via compressive spatiotemporal population receptive fields. Nat Commun 15, 6885 (2024). https://doi.org/10.1038/s41467-024-51243-7

## Data
Data are stored on the Open Science Foundation URL: https://osf.io/rpuhs (main paper) and https://osf.io/e83az/ (supplemental data).

## Dependencies:
SpatiotemporalPRF toolbox: https://github.com/VPNL/spatiotemporalPRFs  
VistaSoft toolbox: https://github.com/vistalab/vistasoft  
NOTE: We are in the process of integrating the spatiotemporal pRF simulation code (stRet toolbox by Insub Kim) into the main vistasoft toolbox. We recommend using Insub Kim's forked repository (https://github.com/KimInsub/vistasoft) in the meantime.

## Overview of analysis pipeline
* simseqRootPath.m
* downloadDataFromOSF.m
* s_simseqAnalysis_overview
    - Preprocess nifti data to mrVista gray
    - Prepare data for modelfits (make ROIs, concatenate condition parfiles, detrend, average, combine, and chop measured time series).
    - Make model predictions given pRFs
    - Fit model predictions to data
    - Find best fitting R^2 for CST pRF exponent (grid fit)
    - General visualization
* s_makeAllManuscriptFigures.m : Recreate manuscript figures 2,3,4,6,7,8 and supplementary figures with data and simulations
* stimulus/
    - folder with code for stimulus generation and running MRI experiment (stim_mri)
* external/
    - folder with borrowed functions from other toolboxes
* analysis/
    - folder with analysis scripts and subfunctions, and figure plotting functions

NB. Participants S1 through S10 in the paper correspond to participants with different naming convention ("subjXX") from a larger dataset shared between multiple projects: [S1: subj01, S2: subj02, S3: subj03, S4: subj07, S5: subj08, S6: subj09, S7: subj10, S8: subj11, S9: subj12, S10: subj13]
