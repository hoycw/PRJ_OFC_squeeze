% Data Selection
an.ROI         = {'all'};             % Channel to be analyzed
an.event_type  = 'S';           % event around which to cut trials
an.trial_lim_s = [-1 2.01];       % window in SEC for cutting trials
an.demean_yn   = 'no';
an.bsln_evnt   = 'A';       % 'A' means baseline correct to entire trial ("all")
an.bsln_type   = 'mad';
an.bsln_lim    = [-0.8 1];  % For 'A', onset is relative to S and offset is relative to D
an.bsln_boots  = 0;

% TFR Parameters
cfg_tfr = [];
cfg_tfr.method     = 'mtmconvol';
cfg_tfr.output     = 'pow';
cfg_tfr.taper      = 'hanning';
cfg_tfr.foi        = [2:40];
cfg_tfr.t_ftimwin  = 4./cfg_tfr.foi;
cfg_tfr.toi        = 'all'; %-0.2:0.004:1.0;
cfg_tfr.keeptrials = 'yes'; % need trials for stats, can average later
cfg_tfr.keeptapers = 'no';
cfg_tfr.pad        = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype    = 'zero';

% Window Logic Check: bsln_lim is within trial_lim_s
if an.bsln_lim(1) < an.trial_lim_s(1) || an.bsln_lim(2) > an.trial_lim_s(2)
    error('an.bsln_lim is outside an.trial_lim_s!');
end
