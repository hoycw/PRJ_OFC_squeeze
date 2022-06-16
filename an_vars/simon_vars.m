%% Simon TFR parameters
an.ROI         = {'all'};             % Channel to be analyzed
an.event_type  = 'S';           % event around which to cut trials
an.trial_lim_s = [-0.25 2.01];       % window in SEC for cutting trials
an.demean_yn   = 'no';
an.bsln_evnt   = 'S';
an.bsln_type   = 'zscore';
an.bsln_lim    = [-0.25 -0.05];
an.bsln_boots  = 0;

cfg=[];
cfg.trials          = 'all';
cfg.keeptrials      = 'yes';
cfg.output          = 'pow';
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.foi             = 2:80;
% data.sampleinfo(1,2)
cfg.toi             = -1:0.02:3;
cfg.pad             = 'maxperlen';
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;
% Frequecy dependent wavelet length %
cfg.t_ftimwin       = 4./cfg.foi;

cfg.baseline     = [-1 0];
cfg.baselinetype = an.bsln_type;    % 'absolute', 'relative' ratio, 'relchange' %, 'normchange', 'db', 'zscore'.
cfg.parameter    = 'powspctrm';

%% more older?
%     databl=data;
%     % Swap over
%     databl.trial=databl.trialbl;
%     databl=rmfield(databl,'trialbl');
%     
%     cfg=[];
%     cfg.trials          = 'all';
%     cfg.keeptrials      = 'yes';
%     cfg.output          = 'pow';
%     cfg.method          = 'mtmconvol';
%     cfg.taper           = 'hanning';
%     cfg.foi             =2:80;
%     % data.sampleinfo(1,2)
%     cfg.toi             = -3:0.02:3;
%     cfg.pad             ='maxperlen'
%     % cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;
%     % Frequecy dependent wavelet length %
%     cfg.t_ftimwin       =4 ./cfg.foi;
%     freqoutblrp             = ft_freqanalysis(cfg, databl);

    % Manual baseline correction due to the non static baseline when event locking %
    % Make copy
%     time=freqoutblrp.time;
%     st=-1;
%     ed=0;
%     [mns Ixst]=min(abs(time-st));
%     [mne Ixed]=min(abs(time-ed));
%     
%     pwrbl=freqoutblrp.powspctrm;
%     pwrblm=mean(pwrbl(:,:,:,Ixst:Ixed),4);
%     pwrlbmr=repmat(pwrblm,1,1,1,size(freqout.powspctrm,4));
    
    % Baseline correction % Find the baseline %
    % Note the time axis of the event locked data is different from the
    % baseline data as I wanted to go back further for the stim locked data %.
    
%     pwr=freqout.powspctrm;
%     pwr=((pwr-pwrlbmr)./pwrlbmr).*100;
    
%     freqoutbl=freqout;
%     freqoutbl.powspctrm=pwr;
    
