% st.stat_lab  = 'thetapk';
% st.find_peak = 1;           % 0/1 for finding patient specific peak

% st.psd_lim   = [0.25 1.5];
% st.peak_sign = 1;
% st.cfreq_lim = [2 7];
% st.bandwidth = 4;

st.stat_evnt = 'S';
st.stat_lim  = [0 2];
st.time_resolved = 1;

% Sliding Window Parameters
st.win_len  = 0.1;
st.win_step = 0.05;
st.win_center = [st.stat_lim(1)+st.win_len/2:st.win_step:st.stat_lim(2)-st.win_len/2];

st.norm_bhv_pred = 'zscore';    % see fn_normalize_predictor
st.norm_nrl_pred = 'zscore';

st.outlier_thresh = 4;          % threshold for tossing outliers in DV (y variable)
st.outlier_stat_id = 'S5t15_bhvz_nrlz_out4';          % threshold for tossing outliers in DV (y variable)

st.use_simon_tfr = 0;
st.toss_same_trials = 1;

