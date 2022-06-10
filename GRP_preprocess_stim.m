%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
clc
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

evnt_lab = 'stim';      % event to time-lock segmented data {'stim' | 'dec'}
trl_lim = [-3 6];       % cut trial data (in sec) relative to event
new_srate = 1000;
trl_lim_samp = trl_lim.*new_srate;
n_choice_trl = 75;

outlier_std_thresh = 3;

%% Load parameters %%
% Columns in SqueezeSubjectSyncSummary.xlsx:
%   Patient, Block1 Neural Data, Block1 Behavioral Data, Blck1 Sync1, Blck1
%   Sync 2, Block2 Neural Data, Block1 Behavioral Data (likely typo?),
%   Blck1 Sync 1, Blck1 Sync 2, SampleRate
DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';

[numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);
SyncDetails = numbers;

%% Process data
for s=1:3
%     clearvars -except s SBJs sbj_pfc_roi evnt_lab trl_lim new_srate trl_lim_samp n_choice_trl Flnum SyncDetails FileDetails
%     close all;
    
    for b_ix = 1:2
        % Set SBJ-specific filenames and sync details
        if b_ix==1
            nrl_fname = fullfile(cd,'DataFiles',FileDetails{s,2});    % Block1 Neural Data
            bhv_fname = fullfile(cd,'DataFiles',FileDetails{s,3});    % Block1 Behavioral Data
            PCS_start    = SyncDetails(s,1);  % Blck1 Sync 1
            biopac_start = SyncDetails(s,2);  % Blck1 Sync 2
        elseif b_ix==2
            nrl_fname = fullfile(cd,'DataFiles',FileDetails{s,6});    % Block2 Neural Data
            bhv_fname = fullfile(cd,'DataFiles',FileDetails{s,7}); % Block2? behavioral data
            PCS_start    = SyncDetails(s,5);    % Blck2? Sync 1
            biopac_start = SyncDetails(s,6);  % Blck2? Sync 2
        else
            error('only 2 blocks!');
        end
        
        % Loading neural data
        load(nrl_fname);
        
        % Synchronise data- CWH: resample to 1 kHz
        % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
        % Biopacs
        srate = SyncDetails(s,7);%422;  % Sample Rate
        n_samp = size(signal,2);
        srorig = srate;
        len_s  = n_samp./srorig;
        dt = 1./srate;
        orig_time_ax = [0:dt:len_s]; orig_time_ax(end)=[];
        nrl_resamp = [];
        for ch_ix = 1:2
            nrl_resamp(ch_ix,:) = resample(signal(ch_ix,:),orig_time_ax,new_srate);
        end
        
        %% Plot PSDs of the data
        [pxx1,f] = pwelch(nrl_resamp(1,:),new_srate,0,[1:100],new_srate);
        [pxx2,f] = pwelch(nrl_resamp(2,:),new_srate,0,[1:100],new_srate);
        figure; loglog(f,pxx1); hold on; loglog(f,pxx2); box off
        xticks([1 4 8 12 20 30:10:100]);
        legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
        title([SBJs{s} ' block ' num2str(b_ix) ' PSDs']);
        
        %% Find and clean up the PCS data
        %     figure; plot(dtrs1(1,:));
        % PCSstr=input('Input timing of the final pulse   ');
        nrl_trim = nrl_resamp(:,PCS_start:end);
        
        
        load(bhv_fname);
        
        % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
        matlab_bsln_time = syncR1.strtm;
        % Find and sync the Biopacs / Matlab data.
        figure; plot(syncR1.data(:,3)); title([SBJs{s} ' block ' num2str(b_ix) ' Biopac data']);
        % Biopstr=input('Input timing of the final pulse');
        
        bhv = result.data;
        
        %% Create Fieldtrip structure
        blk_data{b_ix} = struct;
        % data.fsample=new_srate;
        blk_data{b_ix}.label{1,1}='LFP';
        blk_data{b_ix}.label{2,1}='OFC';
        
        ad1 = biopac_start./500;
        % Have checked the timings of this and it works now.
        % !!! CWH: unclear what these timestamps are, but Simon says above it's okay...
        PCS_comb_bhv_start = matlab_bsln_time+ad1;
        
        % Find where a button was pressed %
        trial_type = zeros([n_choice_trl 1]);
        for trl_ix = 1:n_choice_trl
            if ~isempty(bhv(trl_ix).key)
                trial_type(trl_ix)=bhv(trl_ix).key==37 | bhv(trl_ix).key==39;
            end
        end
        resp_ix = find(trial_type==1);
        
        blk_effort{b_ix}   = nan(size(trial_type));
        blk_stake{b_ix}    = nan(size(trial_type));
        blk_decision{b_ix} = nan(size(trial_type));
        trial_onset        = nan(size(trial_type));
        blk_data{b_ix}.time  = cell(size(trial_type));
        blk_data{b_ix}.trial = cell(size(trial_type));
        for t = 1:length(resp_ix)
            trl_ix = resp_ix(t);
            blk_effort{b_ix}(trl_ix,1)   = bhv(trl_ix).effortIx;
            blk_stake{b_ix}(trl_ix,1)    = bhv(trl_ix).stakeIx;
            trial_onset(trl_ix,1) = bhv(trl_ix).startStim-PCS_comb_bhv_start;  %1.5 seconds until motor mapping revealed.
            
            if bhv(trl_ix).Yestrial ==1
                choice_onset(trl_ix,1) = bhv(trl_ix).YesChoice-PCS_comb_bhv_start;
            else
                choice_onset(trl_ix,1) = bhv(trl_ix).NoChoice-PCS_comb_bhv_start;
            end
            
            blk_decision{b_ix}(trl_ix,1) = bhv(trl_ix).Yestrial;
            
            if strcmp(evnt_lab,'stim')
                start_samp = round(trial_onset(trl_ix)*new_srate) + trl_lim_samp(1);
            elseif strcmp(evnt_lab,'dec')
                start_samp = round(choicestr(trl_ix)*new_srate) + trl_lim_samp(1);
            else
                error('not ready for event besides stim or dec');
            end
            end_samp   = start_samp + trl_lim_samp(2);                         % to end plus 6 seconds after
            blk_data{b_ix}.trial{trl_ix} = nrl_trim(:,start_samp:end_samp);
            
            % Remember this is dependent on sample rate and should be in seconds.
            blk_data{b_ix}.time{trl_ix} = ([start_samp:end_samp]-start_samp+trl_lim_samp(1))./new_srate;
            blk_data{b_ix}.sampleinfo(trl_ix,:) = [start_samp end_samp];
            %     blk_data.trialinfo(trl_ix,1)=trls1(trl_ix).effortIx;
            %     blk_data.trialinfo(trl_ix,2)=trls1(trl_ix).stakeIx;
        end
        
        % Create previous trial regressors
        blk_effort_prv{b_ix} = blk_effort{b_ix};
        blk_effort_prv{b_ix}(end) = [];                     % remove last trial
        blk_effort_prv{b_ix} = [nan; blk_effort_prv{b_ix}]; % add nan to start, shifting by 1
        
        blk_stake_prv{b_ix} = blk_stake{b_ix};
        blk_stake_prv{b_ix}(end) = [];
        blk_stake_prv{b_ix} = [nan; blk_stake_prv{b_ix}];
        
        blk_decision_prv{b_ix} = blk_decision{b_ix};
        blk_decision_prv{b_ix}(end) = [];
        blk_decision_prv{b_ix} = [nan; blk_decision_prv{b_ix}];
        
        clear syncR1 syncR2 result signal nrl_resamp nrl_trim orig_time_ax power
%         close all
    end
    
    %% Combine data across blocks
    % Combine behavioral data
    effort   = [blk_effort{1}; blk_effort{2}];
    stake    = [blk_stake{1}; blk_stake{2}];
    decision = [blk_decision{1}; blk_decision{2}];
    effort_prv   = [blk_effort_prv{1}; blk_effort_prv{2}];
    stake_prv    = [blk_stake_prv{1}; blk_stake_prv{2}];
    decision_prv = [blk_decision_prv{1}; blk_decision_prv{2}];
    
    % Combine neural data
    data = blk_data{1};
    for trl_ix=1:n_choice_trl
        data.trial{trl_ix+n_choice_trl} = blk_data{2}.trial{trl_ix};
        data.time{trl_ix+n_choice_trl}  = blk_data{2}.time{trl_ix};
%         data.trialbl{trl_ix+n_choice_trl} = blk_data{2}.trialbl{n_choice_trl};
    end
    
    %% Find where the was a missing trial (not completed).
    empty_ix = find(effort==0);
    if ~isempty(empty_ix); error('why are there still empty trials after excluding above?'); end
    
%     effort(whrempty)=[];
%     stake(whrempty)=[];
%     decision(whrempty)=[];
%     data.trial(whrempty)=[];
%     data.time(whrempty)=[];
%     data.trialbl(whrempty)=[];
    
    % data.sampleinfo(76:150,:)=data2.sampleinfo;
    
    %% Find and reject outliers
    % Get standard deviation and normalize per channel
    trial_std = nan([size(data.trial,1) length(data.label)]);
    for trl_ix = 1:size(data.trial,1)
        tmp = data.trial{trl_ix};
        trial_std(trl_ix,:)  = std(tmp,1,2);
%         trialmean(trl_ix,:) = mean(abs(tmp),2);    % why are you taking abs here?
    end
    trial_std_norm = trial_std./mean(trial_std);
    
    % figure;
    % subplot(2,1,1); plot(trialstdN(:,1));
    % subplot(2,1,2); plot(trialmean);
    
    % Find trials with outlier standard deviations
    bad_trl_ix = [];
    for ch_ix = 1:2
        bad_trl_ix = [bad_trl_ix; find(trial_std_norm(:,ch_ix)>outlier_std_thresh)];
    end
    bad_trl_ix = unique(bad_trl_ix);
    if isempty(bad_trl_ix)
        fprintf('%s: No bad trials identified.\n',SBJs{s});
    else
        fprintf(2,'%s: Removing %d trials!\n',SBJs{s},length(bad_trl_ix));
    end
    
    % Remove outlier trials
    data.trial(bad_trl_ix) = [];
%     data.trialbl(bad_trl_ix) = [];
    data.time(bad_trl_ix) = [];
    data.sampleinfo(bad_trl_ix,:) = [];
    
    effort(bad_trl_ix)   = [];
    stake(bad_trl_ix)    = [];
    decision(bad_trl_ix) = [];
    effort_prv(bad_trl_ix)   = [];
    stake_prv(bad_trl_ix)    = [];
    decision_prv(bad_trl_ix) = [];
    
    %% BEHAVIORAL MODELLING %%
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    effortlevels = [0.16 0.32 0.48 0.64 0.80];
    stakelevels  = [1 4 7 10 13];
    
    effortfr=zeros(size(effort,1),1);
    for g=1:5
        wh=find(effort==g);
        effortfr(wh)=effortlevels(g);
    end
    
    for g=1:5;
        stake2(stake==g)=stakelevels(g);
    end
    
    stake=stake2;
    
    % Fit the behavior. Minimise the difference between the probability and the decision
    decisionfun=@(p) norm( (exp(p(1)*(stake-(p(2)*(effortfr).^2))) ./ (exp(p(1)) + exp(p(1)*(stake-(p(2)*(effortfr).^2))))) - decision);
    [par fit]=fminsearch(decisionfun, [1,1]);
    
    SV=@(k) stake-(k*(effortfr).^2);
    EFF=@(k) (k*(effortfr).^2);
    SubjVal1=SV(par(2));
    EFFs=EFF(par(2));
           
    figure;
    subplot(3,1,1);
    scatter(stake,SV(par(2)))
    xlabel('Apples'); ylabel('SV');
    subplot(3,1,2);
    scatter(effortfr,SV(par(2)))
    xlabel('Effort - Proportion max'); ylabel('SV');
    subplot(3,1,3);
    scatter(SV(par(2)),(exp(par(1)*(stake-(par(2)*(effortfr).^2))) ./ (exp(par(1)) + exp(par(1)*(stake-(par(2)*(effortfr).^2))))))
    xlabel('SV'); ylabel('Probabilty of acceptance');

    % Plot a regression line %
    SVtmp=SV(par(2));
    ProbAccept=(exp(par(1)*(stake-(par(2)*(effortfr).^2))) ./ (exp(par(1)) + exp(par(1)*(stake-(par(2)*(effortfr).^2)))))
    
%     figure; scatter(SVtmp,ProbAccept,'.')  

%     sigparam=sigm_fit(SVtmp,ProbAccept,1)
    
%     [p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
%     PotValues=[-10:0.1:10];
%     Pout=polyval(p,PotValues);
    
%     figure; plot(PotValues,Pout)
    
    %% Test behavioral modelling.
    
    figure;
    subplot(2,1,1);
    scatter(stake, SubjVal1);
    title('Reward')
    subplot(2,1,2);
    scatter(effortfr, SubjVal1);
    title('Effort');
    
    % Regression
    % Create regression
    RG=[ones(size(stake)),stake,effortfr];
    betas=regress(SubjVal1,RG)%% FIELDTRIP TF ANALYSIS
    
    % Start by ignoring the autonomic data and just do a correlation of the
    % neural features by the decision features.
    
    
    
    
    %% Time Frequency analysis on the dataset.
    cfg=[];
    cfg.trials          = 'all';
    cfg.keeptrials      = 'yes';
    cfg.output          = 'pow';
    cfg.method          = 'mtmconvol';
    cfg.taper           = 'hanning';
    cfg.foi             =2:80;
    % data.sampleinfo(1,2)
    cfg.toi             = -3:0.02:3;
    cfg.pad             ='maxperlen'
    % cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;
    % Frequecy dependent wavelet length %
    cfg.t_ftimwin       =4 ./cfg.foi;
    freqout             = ft_freqanalysis(cfg, data);
    
    cfg.baseline     = [-1 0]
    cfg.baselinetype = 'zscore'    % 'absolute', 'relative' ratio, 'relchange' %, 'normchange', 'db', 'zscore'.
    cfg.parameter    = 'powspctrm';
    freqoutbl=ft_freqbaseline(cfg,freqout)

    
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
    
    %% SAVE data
    AllData=struct;
    AllData.desp=strcat(FileDetails{s,1},'Stimulus_Locked'); % Patient code
    AllData.raw=data;
    AllData.TF=freqout;
    AllData.TFbl=freqoutbl;
    AllData.TFcfg=cfg;
    AllData.exp.stake=stake;
    AllData.exp.effort=effortfr;
    AllData.exp.EFFs=EFFs;
    AllData.exp.decision=decision;
    AllData.exp.SV=SubjVal1;
    AllData.beh.par=par;
    AllData.beh.fit=fit;
    AllData.expS.stakeS=stakeS;
    AllData.expS.decisionS=decisionS;
    AllData.expS.effortS=effort_prv;
    
    outfold='C:\Users\slittle\Box\Research\Motivation\OFC_recordings_squeeze\Analysis\Group\Preprocess\Output_files_shifted_behavior';
    
%     save(fullfile(outfold,strcat(FileDetails{Flnum,1},'Stimulus_Locked')),'AllData');
    
%     cfg = [];
%     figure;
%     subplot(2,1,1);
%     cfg.channel      = 'OFC';
%     ft_singleplotTFR(cfg,freqout);
%     xlim([-1.5 1.5]);
%     subplot(2,1,2);
%     cfg.channel      = 'LFP';
%     ft_singleplotTFR(cfg,freqout);
%     xlim([-1.5 1.5]);
    
    cfg = [];
    figure;
    subplot(2,1,1);
    cfg.channel      = 'OFC';
    ft_singleplotTFR(cfg,freqoutbl);
    xlim([-0.5 2]);
    ylim([ 2 35])
    subplot(2,1,2);
    cfg.channel      = 'LFP';
    ft_singleplotTFR(cfg,freqoutbl);
    xlim([-0.5 2]);
    ylim([ 2 35])
    
    %% Statistical testing %%
    
    % Start by examining for a change in the power with the event %
    
    % ft_statistics_montecarlo
    
    % Start by examining the reward and effort separately %
    % compute statistics with ft_statfun_indepsamplesregrT
    
    cfg = [];
    cfg.channel          = 'OFC';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    % cfg.statistic        = 'ft_statfun_depsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    % This should be 1 - 3 as this actually include the post presentation
    % period only.
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = stake';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statOFCstake = ft_freqstatistics(cfg, freqoutbl);
    
    cfg = [];
    cfg.channel          = 'OFC';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    % cfg.statistic        = 'ft_statfun_depsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = EFF(par(2))';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statOFCeffort = ft_freqstatistics(cfg, freqoutbl);
    
    cfg = [];
    cfg.channel          = 'OFC';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    % cfg.statistic        = 'ft_statfun_depsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = SubjVal1';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statOFCSV = ft_freqstatistics(cfg, freqoutbl);
    
    cfg = [];
    cfg.channel          = 'LFP';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = stake';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statLFPstake = ft_freqstatistics(cfg, freqoutbl);
    
    cfg = [];
    cfg.channel          = 'LFP';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = EFF(par(2))';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statLFPeffort = ft_freqstatistics(cfg, freqoutbl);
    
    cfg = [];
    cfg.channel          = 'LFP';
    cfg.statistic        = 'ft_statfun_indepsamplesregrT';
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.numrandomization = 1000;
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.correcttail      ='alpha';
    cfg.frequency        = [2 30];
    cfg.latency          = [0 1.5];
    
    n1 = size(stake,1);
    design(1,1:n1)       = SubjVal1';
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    statLFPSV = ft_freqstatistics(cfg, freqoutbl);
    
    % PLOT STATS RESULTS
    
    figure;
    subplot(2,3,1);
    cfg           = [];
    cfg.channel   = {'OFC'};
    cfg.parameter = 'stat';
    cfg.colormap  = parula;
    cfg.ylim      = [2 100];
    cfg.xlim      = [-1.5 1.5];
    % cfg.zlim      = [ ] ;
    cfg.marker    ='off';
    cfg.style     = 'fill';
    cfg.comment   = 'off';
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    cfg.colorbar  = 'yes';
    ft_singleplotTFR(cfg,statOFCstake);
    
    subplot(2,3,2);
    cfg.channel   = {'OFC'};
    ft_singleplotTFR(cfg,statOFCeffort);
    
    subplot(2,3,3);
    cfg.channel   = {'OFC'};
    ft_singleplotTFR(cfg,statOFCSV);
    
    subplot(2,3,4);
    cfg.channel   = {'LFP'};
    ft_singleplotTFR(cfg,statLFPstake);
    
    subplot(2,3,5);
    cfg.channel   = {'LFP'};
    ft_singleplotTFR(cfg,statLFPeffort);
    
    subplot(2,3,6);
    cfg.channel   = {'LFP'};
    ft_singleplotTFR(cfg,statLFPSV);
end