# PFC squeeze paper code repository
This code repository is for all preprocessing, analysis, and plotting scripts related to Hoy, de Hemptinne, et al. "Beta and theta oscillations track effort and previous reward in human basal ganglia and prefrontal cortex during decision making" (biorxiv 2024, in review)

## Data and Experimental Paradigm:
Experimental paradigm as in Saleh et al., (2021; Brain;  doi:10.1093/brain/awab013). 
All raw and preprocessed datasets are provided under Creative Commons Attrubition 4.0 license on OSF at DOI 10.17605/OSF.IO/JQ6Z4.
This repo is not citable (yet).

## Manuscript: 
Manuscript is in review.
Initial biorxiv preprint: [https://www.biorxiv.org/content/10.1101/2022.12.07.519496v1]. 
Code written by Colin W. Hoy and Simon Little.

## Dependencies
OS: MacBook Pro running OS 11.6.1; MATLAB version R2023b
  - External toolboxes:
    - Fieldtrip (version 7ba023993): <http://www.fieldtriptoolbox.org/>

## Overview of Code
Core scripts to reproduce analyses in the manuscript should be run in order. Note that "orig_srate" and "osr" are the final versions of the code reported in the manuscript. The biorxiv preprint and original submission used an upsampling preprocessing step for the neural data to match the Biopac sampling rate of the behavioral data stream, but upon reviewer request, we removed this step, which improved the strength of the (unchanged) results. Thus, "osr" versions are reported as final.
1. GRP00 takes raw data, aligns, extracts trial-wise neural and behavioral data, tosses bad trials, saves
2. GRP01a- TFR filtering, baseline correction, save
3. GRP01b- plots and save TFRs
4. GRP02- plot TFRs and find subject-specific peaks (e.g., low beta)
5. GRP03- average trial-wise power and combine in a table with behavioral predictors for LMMs
6. GRP04- run LMM models and plot results for different power variables (theta, beta, connectivity)
7. GRP05- run LMM models for RT data
8. GRP06- plot behavior (this is the one that is somewhat out of order)
9. GRP0304- combines GRP03 and GRP04 to run the LMMs at each time point, purely to visualize the evolution of the effects
10. PFC05- run LMMs for PFC stim behavior

## Overview of repo structure
    - SBJ_vars.m contains information about each subject for consistent processing
    - an_vars/ contains parameters for signal processing pipelines. Final verison for the paper is an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'
    - stat_vars/ contains parameters for statistical analyses. Final version for the paper is stat_id = 'S5t15_bhvz_nrl0_out3'
    - plt_vars/ contains plotting parameters, which vary depending on the figure. plt_id = 'ts_S2t2_evnts_sigLine' is commonly used.

## Supplemental Analyses
- GRP03b- demonstrates the effect of thresholds for rejecting outlier trials in single-trial power
- "perm" versions of GRP04 and GRP0304 used non-parametric permutation statistics to verify our main findings (see manuscript Methods).
- GRP06b has additional attempts at integrating previous reward into the computational model of behavior. These largely failed, as novel model development is best conducted on large datasets that allow proper testing and validation.
- Additional plotting scripts to understand the data:
  - GRP01c scripts show power differences for various splits of the behavioral variables (reviewer figures).
  - GRP01d was used to assess reliability of evoked potentials.
  - GRP01e was used to check for high frequency artifacts.
- Preliminary connectivity analyses were run, but these results were noisy and deemed unreliable. There was a hint of PFC-BG theta amplitude correlations tracking previous reward, but again, we decided this was not worth reporting.
  - GRP07a- compute PFC-BG connectivity
  - GRP07b- plot connectivity
  - GRP08- bandpass data (not using GRP07 output!) and average trial-wise connectivity, build table with predictors
  - GRP09- run LMM models for connectivity
  - GRP10- circular-linear correlation between phase and behavioral predictors

