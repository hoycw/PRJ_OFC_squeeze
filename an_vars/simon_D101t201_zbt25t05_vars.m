%% Simon TFR parameters
an.ROI         = {'all'};             % Channel to be analyzed
an.event_type  = 'D';           % event around which to cut trials
an.trial_lim_s = [-1.01 2.01];       % window in SEC for cutting trials
an.demean_yn   = 'no';
an.bsln_evnt   = 'S';
an.bsln_type   = 'zboot';
an.bsln_lim    = [-0.25 -0.05];
an.bsln_boots  = 500;

cfg_tfr=[];
cfg_tfr.trials          = 'all';
cfg_tfr.keeptrials      = 'yes';
cfg_tfr.output          = 'pow';
cfg_tfr.method          = 'mtmconvol';
cfg_tfr.taper           = 'hanning';
cfg_tfr.foi             = [2:0.5:8 8:20 22:2:40];
% cfg_tfr.foi             = 2:80;
% data.sampleinfo(1,2)
cfg_tfr.toi             = 'all';%-3:0.02:3;
cfg_tfr.pad             = 'maxperlen';
% cfg_tfr.t_ftimwin    = ones(length(cfg.foi),1).*1;
% Frequecy dependent wavelet length %
cfg_tfr.t_ftimwin       = 4./cfg_tfr.foi;

% cfg.baseline     = an.bsln_lim;
% cfg.baselinetype = an.bsln_type;    % 'absolute', 'relative' ratio, 'relchange' %, 'normchange', 'db', 'zscore'.
% cfg.parameter    = 'powspctrm';

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
    
