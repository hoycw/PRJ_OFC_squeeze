% Data Selection
an.ROI         = {'all'};             % Channel to be analyzed
an.event_type  = 'S';           % event around which to cut trials
an.trial_lim_s = [-1 2.01];       % window in SEC for cutting trials
an.demean_yn   = 'no';
an.bsln_evnt   = 'S';
an.bsln_type   = 'mad';
an.bsln_lim    = [-0.8 0];
an.bsln_boots  = 0;

% TFR Parameters
cfg_tfr = [];
cfg_tfr.method     = 'mtmconvol';
cfg_tfr.output     = 'pow';
cfg_tfr.taper      = 'hanning';
cfg_tfr.foi        = [80:120];
cfg_tfr.t_ftimwin  = 4./cfg_tfr.foi;
cfg_tfr.toi        = an.trial_lim_s(1):0.004:an.trial_lim_s(2);
cfg_tfr.keeptrials = 'yes'; % need trials for stats, can average later
cfg_tfr.keeptapers = 'no';
cfg_tfr.pad        = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype    = 'zero';

% Use preprocessed data with original sampling rate
an.orig_srate = 1;

% Window Logic Check: bsln_lim is within trial_lim_s
if an.bsln_lim(1) < an.trial_lim_s(1) || an.bsln_lim(2) > an.trial_lim_s(2)
    error('an.bsln_lim is outside an.trial_lim_s!');
end
