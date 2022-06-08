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
    clearvars -except s SBJs sbj_pfc_roi Flnum SyncDetails FileDetails
    close all;
    
    % Loading neural data
    fname=fullfile(cd,'DataFiles',FileDetails{s,2});    % Block1 Neural Data
    load(fname);
    
    % Synchronise data- CWH: resample to 1 kHz
    % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
    % Biopacs
    srate = SyncDetails(s,7);%422;  % Sample Rate
    new_srate = 1000;
    n_samp = size(signal,2);
    srorig = srate;
    len_s  = n_samp./srorig;
    dt = 1./srate;
    orig_time_ax=[0:dt:len_s]; orig_time_ax(end)=[];
    dtrs1=[];
    for ch_ix = 1:2
        dtrs1(ch_ix,:)=resample(signal(ch_ix,:),orig_time_ax,new_srate);
    end
    
    %% Plot PSDs of the data
    [pxx1,f] = pwelch(dtrs1(1,:),mFs,0,[1:100],mFs);
    [pxx2,f] = pwelch(dtrs1(2,:),mFs,0,[1:100],mFs);
    figure; loglog(f,pxx1); hold on; loglog(f,pxx2); box off
    legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
    
    %% Find and clean up the PCS data
%     figure; plot(dtrs1(1,:));
    % PCSstr=input('Input timing of the final pulse   ');
      
    PCSstr = SyncDetails(s,1);  % Blck1 Sync 1
    
    dtrst1=dtrs1(:,PCSstr:end);
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,3});    % Block1 Behavioral Data
    load(fname);

    mlbltm1=syncR1.strtm;
    % Find and sync the Biopacs / Matlab data.
    figure; plot(syncR1.data(:,3));
    % Biopstr=input('Input timing of the final pulse');
    Biopstr1=SyncDetails(s,2);  % Blck1 Sync 2
    
    trls1=result.data;
    
    %% Create Fieldtrip structure
    data=struct;
    % data.fsample=mFs;
    data.label{1,1}='LFP';
    data.label{2,1}='OFC';
    
    ad1=Biopstr1./500;
    % Have checked the timings of this and it works now.
    PCSst1=mlbltm1+ad1;
    
    % Find where a button was pressed %
    trialtype=[];
    for trl_ix = 1:75 % !!! why only to 75? there are 85 trials in trls1
        if ~isempty(trls1(trl_ix).key)
        trialtype(trl_ix)=trls1(trl_ix).key==37 | trls1(trl_ix).key==39;
        end
    end
    whrtrl=find(trialtype==1);
    
    for trl_ix=whrtrl
        effort(trl_ix,1)=trls1(trl_ix).effortIx;
        stake(trl_ix,1)=trls1(trl_ix).stakeIx;
        trialstr(trl_ix,1)=trls1(trl_ix).startStim-PCSst1;  %1.5 seconds until motor mapping revealed.
        
%         if trls1(trl_ix).Yestrial ==1
%             choicestr(trl_ix,1)=trls1(trl_ix).YesChoice-PCSst1;
%         else
%             choicestr(trl_ix,1)=trls1(trl_ix).NoChoice-PCSst1;
%         end
%         
        decision(trl_ix,1)=trls1(trl_ix).Yestrial;
        
        st=round(trialstr(trl_ix)*1000)-3000; % grab start plus 3 seconds before
        ed=st+6000;                         % to end plus 6 seconds after
        data.trial{trl_ix}=dtrst1(:,st:ed);
        
%         st=round(choicestr(trl_ix)*1000)-3000;
%         ed=st+6000;
%         data.trial{trl_ix}=dtrst1(:,st:ed);
        
        data.trial{trl_ix}=dtrst1(:,st:ed);
        % Remember this is dependent on sample rate and should be in seconds.
        data.time{trl_ix}=([st:ed]-st-3000)./1000;
        data.sampleinfo(trl_ix,:)=[st ed];
        %     data.trialinfo(trl_ix,1)=trls1(trl_ix).effortIx;
        %     data.trialinfo(trl_ix,2)=trls1(trl_ix).stakeIx;
        
        
    end
    
    %% Load Block 2 data
    clear syncR1 syncR2 result signal
    close all
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,6}); % Block2 Neural Data
    load(fname);
       
    % Synchronise data
    % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
    % Biopacs
    n_samp = size(signal,2);
    srorig = srate;
    len_s  = n_samp./srorig;
    dt = 1./srate;
    orig_time_ax=[0:dt:len_s]; orig_time_ax(end)=[];
    dtrs2=[];
    for ch_ix = 1:2
        dtrs2(ch_ix,:)=resample(signal(ch_ix,:),orig_time_ax,new_srate);
    end
    
    % Plot PSDs of the data
    [pxx1,f] = pwelch(dtrs2(1,:),mFs,0,[1:100],mFs);
    [pxx2,f] = pwelch(dtrs2(2,:),mFs,0,[1:100],mFs);
    figure; loglog(f,pxx1); hold on; loglog(f,pxx2); box off
    legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
    
    % Find and clean up the PCS data
%     figure; plot(dtrs2(1,:));
    % PCSstr=input('Input timing of the final pulse   ');
    % Note that for this subject - using the first pulse as the end is missing. (Ignore the first
    % deflection which is 300ms before the first pulse, this seems to be a
    % turing on artefact and is not visible in the EMG pulse')
    
    PCSstr=SyncDetails(s,5);    % Blck2? Sync 1
    
    dtrst2=dtrs2(:,PCSstr:end);
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,7}); % Block2? behavioral data
    load(fname);
    
    mlbltm2=syncR1.strtm;
    % Find and sync the Biopacs / Matlab data.
    figure; plot(syncR1.data(:,3));
    % Biopstr=input('Input timing of the final pulse');
    % Here using the first pulse but be careful to ignore the ECG artefact %
    
    Biopstr2=SyncDetails(s,6);  % Blck2? Sync 2
    
    trls2=result.data;
    
    % Create Fieldtrip structure
    data2=struct;
    % data.fsample=mFs;
    data2.label{1,1}='LFP';
    data2.label{2,1}='OFC';
    
    % Need to use the sampling rate here of the Biopacs
    ad2=Biopstr2./500;
    PCSst2=mlbltm2+ad2;
    
    % Find the trials where a key was pressed.
    trialtype=[];
    for trl_ix=1:75
        if ~isempty(trls2(trl_ix).key)
        trialtype(trl_ix)=trls2(trl_ix).key==37 | trls2(trl_ix).key==39;
        end
    end
    whrtrl=find(trialtype==1);
    
    for trl_ix=whrtrl
        effort2(trl_ix,1)=trls2(trl_ix).effortIx;
        stake2(trl_ix,1)=trls2(trl_ix).stakeIx;
        trialstr2(trl_ix,1)=trls2(trl_ix).startStim-PCSst2;  %1.5 seconds until motor mapping revealed.
        
%         if trls2(trl_ix).Yestrial ==1
%             choicestr2(trl_ix,1)=trls2(trl_ix).YesChoice-PCSst2;
%         else
%             choicestr2(trl_ix,1)=trls2(trl_ix).NoChoice-PCSst2;
%         end
        
        decision2(trl_ix,1)=trls2(trl_ix).Yestrial;
        
        st=round(trialstr2(trl_ix)*1000)-3000;
        ed=st+6000;
        data2.trial{trl_ix}=dtrst2(:,st:ed);
        
%         st=round(choicestr2(trl_ix)*1000)-3000;
%         ed=st+6000;
%         data2.trial{trl_ix}=dtrst2(:,st:ed);
        
        % Remember this is dependent on sample rate and should be in seconds.
        data2.time{trl_ix}=([st:ed]-st-3000)./1000;
        data2.sampleinfo(trl_ix,:)=[st ed];
        %     data2.trialinfo(trl_ix,1)=trls2(trl_ix).effortIx;
        %     data2.trialinfo(trl_ix,2)=trls2(trl_ix).stakeIx;
    end
    
    close all;
    clear syncR1 syncR2 result signal
    
    %% Combine together
    
    % Create previous trial regressors
    effortS=effort;
    effortS2=effort2;
    effortS(end)=[];            % remove last trial
    effortS=[nan;effortS];      % add nan to start, shifting by 1
    effortS2(end)=[];
    effortS2=[nan;effortS2];    
    
    stakeS=stake;
    stakeS2=stake2;
    stakeS(end)=[];
    stakeS=[nan;stakeS];
    stakeS2(end)=[];
    stakeS2=[nan;stakeS2];
    
    decisionS=decision;
    decisionS(end)=[];
    decisionS=[nan;decisionS];
    decisionS2=decision2;
    decisionS2(end)=[];
    decisionS2=[nan;decisionS2];

    % Combine both blocks
    effortS=[effortS; effortS2];
    decisionS=[decisionS;decisionS2];
    stakeS=[stakeS;stakeS2];
    
    % Just add data2 to data 55
    effort=[effort;effort2];
    stake=[stake; stake2];
    decision=[decision;decision2];

    for trl_ix=1:75
        data.trial{trl_ix+75}=data2.trial{trl_ix};
        data.time{trl_ix+75}=data2.time{trl_ix};
%         data.trialbl{trl_ix+75}=data2.trialbl{75};
    end
    
    %% Find where the was a missing trial (not completed).
    whrempty=find(effort==0);
    
    effort(whrempty)=[];
    stake(whrempty)=[];
    decision(whrempty)=[];
    data.trial(whrempty)=[];
    data.time(whrempty)=[];
%     data.trialbl(whrempty)=[];
    
    % data.sampleinfo(76:150,:)=data2.sampleinfo;
    
    %% Find and reject outliers
    for trl_ix=1:size(data.trial,2)
        tmp=data.trial{trl_ix};
        trialstd(trl_ix,:)=std(tmp,1,2);
        trialmean(trl_ix,:)=mean(abs(tmp),2);    % why are you taking abs here?
    end
    
    mnstd=mean(trialstd);
    trialstdN=trialstd./mnstd;
    
    % figure;
    % subplot(2,1,1); plot(trialstdN(:,1));
    % subplot(2,1,2); plot(trialmean);
    
    whrc=[];
    for ch_ix = 1:2
        whr=find(trialstdN(:,ch_ix)>3);
        whrc=[whrc;whr];
    end
    
    whrc=unique(whrc);
    
    % Now delete trial whrc on in all data.
    
    data.trial(whrc)=[];
%     data.trialbl(whrc)=[];
    data.time(whrc)=[];
    data.sampleinfo(whrc,:)=[];
    
    effort(whrc)=[];
    stake(whrc)=[];
    decision(whrc)=[];
    
    %% BEHAVIORAL MODELLING %%
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    effortlevels =[0.16 0.32 0.48 0.64 0.80];
    stakelevels  =[1 4 7 10 13];
    
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
    AllData.expS.effortS=effortS;
    
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