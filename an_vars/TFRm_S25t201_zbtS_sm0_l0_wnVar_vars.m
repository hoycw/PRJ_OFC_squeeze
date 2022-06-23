an.evnt_lab    = 'S';           % event around which to cut trials
% trial_lim_s will NOT be full of data! the first and last t_ftimwin/2 epochs will be NaNs
an.trial_lim_s = [-0.25 2.01];      % window in SEC for cutting trials
an.demean_yn   = 'no';             % z-score for HFA instead
an.bsln_evnt   = 'S';
an.bsln_type   = 'zboot';
an.bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
an.bsln_boots  = 1000;              % repetitions for non-parametric stats

% HFA Calculations
an.foi_center  = [2:0.5:8 8:20 22:2:40];
an.octave      = 3/4;              % Frequency resolution
an.foi_min     = 2^(-an.octave/2)*an.foi_center;
an.foi_max     = 2^(an.octave/2)*an.foi_center;
an.foi         = (an.foi_min+an.foi_max)/2;
an.delta_freq  = an.foi_max-an.foi_min;
an.delta_time  = 0.05;
an.n_taper_all = max(1,round(an.delta_freq.*an.delta_time-1));   %number of tapers for each frequency
an.foi_center  = round(an.foi_center*10)/10;          %convert to float?
an.delta_freq_true = (an.n_taper_all+1)./an.delta_time; % total bandwidth around

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'mtmconvol';
cfg_tfr.taper        = 'dpss';
cfg_tfr.tapsmofrq    = an.delta_freq_true./2;                  %ft wants half bandwidth around the foi
cfg_tfr.keeptapers   = 'no';
cfg_tfr.pad          = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = an.foi_center;                          % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_tfr.t_ftimwin    = ones(length(cfg_tfr.foi),1).*an.delta_time;    % length of time window; 0.5 sec, could be n_cycles./foi for n_cylces per win
cfg_tfr.toi          = 'all';%-buff_lim(1):0.1:1.5;         % time window centers
cfg_tfr.keeptrials   = 'yes';                               % must be 'yes' for stats
% cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;

% Log Tranform
an.log_yn = 0;

% Outlier Rejection
% an.outlier_std_lim = 6;

% Cleaning up power time series for plotting
an.smooth_pow_ts = 0;
an.lp_yn       = 'no';
an.lp_freq     = 10;
an.hp_yn       = 'no';
an.hp_freq     = 0.5;

% Resampling
an.resample_ts   = 0;
% an.resample_freq = 250;

