# PRJ_OFC_squeeze

script order:
GRP00 takes raw data, aligns, extracts trial-wise neural and behavioral data, tosses bad trials, saves
GRP01a- TFR filtering, baseline correction, save
GRP01b- plots and save TFRs
GRP02- plot TFRs and find subject-specific peaks
GRP03- average trial-wise power and combine in a table with behavioral predictors for LMMs
GRP04- run LMM models and plot results for different power variables (theta, beta, connectivity)
GRP05- run LMM models for RT data
GRP06- plot behavior (this is the one that is somewhat out of order)
GRP07a- compute PFC-BG connectivity
GRP07b- plot connectivity
GRP08- bandpass data (not using GRP07 output!) and average trial-wise connectivity, build table with predictors
GRP09- run LMM models for connectivity
PFC05- run LMMs for PFC stim behavior
