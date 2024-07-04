# simseqPRF
Code accompanying the paper "Rethinking simultaneous suppression in visual cortex via compressive spatiotemporal population receptive fields" by Kupers, Kim, & Grill-Spector

## Goal
The goal of this project is to operationalize and elucidate how simultaneous and sequential visual stimuli are processed within population receptive fields, in space and time, and generate a lower response for simultaneous over sequential presentations.

## Paper
Title: Rethinking simultaneous suppression in visual cortex via compressive spatiotemporal population receptive fields.
Authors: Kupers, Kim, Grill-Spector
Year: 2024
Journal: Nature Communications
DOI: XXX

## Data
Data are stored on the Open Science Foundation URL: https://osf.io/rpuhs (main paper) and https://osf.io/e83az/ (supplemental data).

## Dependencies:
SpatiotemporalPRF toolbox: https://github.com/VPNL/spatiotemporalPRFs  
VistaSoft toolbox: https://github.com/vistalab/vistasoft  
--> NB current version of VistaSoft might still be in the process of implementing spatiotemporal pRF simulation code. In the meantime, you can use Insub Kim's forked repository here: https://github.com/KimInsub/vistasoft

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
* s_makeAllManuscriptFigures.m : Recreate manuscript figures 2,3,4,6,7,8 and supplementary figures
* stimulus/
    - folder with code for stimulus generation and running MRI experiment (stim_mri)
* external/
    - folder with code from other toolboxes
* analysis/
    - folder with analysis scripts and subfunctions, and figure plotting functions

NB. Participants S1 through S10 in the paper correspond to participants with different naming convention ("subjXX") from a larger dataset shared between multiple projects: [S1: subj01, S2: subj02, S3: subj03, S4: subj07, S5: subj08, S6: subj09, S7: subj10, S8: subj11, S9: subj12, S10: subj13]
