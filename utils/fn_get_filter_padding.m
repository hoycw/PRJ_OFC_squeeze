function [trial_lim_s_pad] = fn_get_filter_padding(cfg_tfr,trial_lim_s,bsln_lim)
%% Determine Filter Padding
%   At a minimum, trial_lim_s must extend 1/2*max(filter_window) prior to 
%   the first data point to be estimated to avoid edges in the filtering window.
%   (ft_freqanalysis will return NaN for partially empty windows, e.g. an edge pre-trial,
%   but ft_preprocessing would return a filtered time series with an edge artifact.)
%   Also note this padding buffer should be at least 3x the slowest cycle
%   of interest.
if strcmp(cfg_tfr.method,'mtmconvol')
    % Cover at least 3 cycles of slowest frequency
    pad_len = 0.5*max(cfg_tfr.t_ftimwin)*3;
elseif strcmp(cfg_tfr.method,'wavelet')
    % add 250 ms as a rule of thumb, or longer if necessary
    pad_len = 0.5*max([cfg_tfr.width/min(cfg_tfr.foi) 0.25]);
else
    error(['Unknown cfg_tfr.method: ' cfg_tfr.method]);
end

% Cut data to bsln_lim to be consistent across S and R/F locked (confirmed below)
%   Add extra 10 ms just because trimming back down to trial_lim_s exactly leaves
%   one NaN on the end (smoothing across that will NaN out everything)
trial_lim_s_pad = [min(bsln_lim)-pad_len trial_lim_s(2)+pad_len+0.01];

% Check that baseline will be included in data cut to trial_lim_s
if trial_lim_s(1) < bsln_lim(1)
    error('ERROR: trial_lim_s does not include bsln_lim!');
end

% Check that trial_lim_s includes full baseline (e.g., zbtA: z-score over all data)
if trial_lim_s_pad(2) < bsln_lim(2)+pad_len+0.01
    trial_lim_s_pad(2) = bsln_lim(2)+pad_len+0.01;
end

