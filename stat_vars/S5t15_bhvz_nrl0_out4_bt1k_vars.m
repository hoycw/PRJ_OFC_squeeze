% st.stat_lab  = 'thetapk';
% st.find_peak = 1;           % 0/1 for finding patient specific peak

% st.psd_lim   = [0.25 1.5];
% st.peak_sign = 1;
% st.cfreq_lim = [2 7];
% st.bandwidth = 4;

st.stat_evnt = 'S';
st.stat_lim  = [0.5 1.5];

st.norm_bhv_pred = 'zscore';    % see fn_normalize_predictor
st.norm_nrl_pred = 'none';

st.outlier_thresh = 4;          % threshold for tossing outliers in DV (y variable)
st.nboots = 1000;

st.use_simon_tfr = 0;
st.toss_same_trials = 1;

