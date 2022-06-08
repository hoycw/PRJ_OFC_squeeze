%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked_PFC01_compatible.m
clear all
close all
clc
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

%% Load data %%

DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';

[numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);
SyncDetails = numbers;

%% Process data
for s=4
    SBJ = FileDetails{s,1};
%     1:size(FileDetails,1)
    clearvars -except s SBJs sbj_pfc_roi SyncDetails FileDetails
    trlNan=[];
    close all;
    fname=fullfile(cd,'DataFiles',FileDetails{s,2});    
    load(fname);
    
    % PFC01 - put it back into signal.
    if strcmp(SBJ,'PFC01') && ~exist('signal','var')
        signal(1,:)=datar(:,1)';        % 1 = +2 -0      % 2 = +1 -3
        signal(2,:)=datar(:,3)';        % 3 = +9 -8    % 4 = +11 -10
    end
    
    % Synchronise data
    % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
    % Biopacs
    Fs=SyncDetails(s,7);
    mFs=1000;
    dtrs1=[];
    for g=1:2
        lgn=size(signal,2);
        srorig=Fs;
        tm=lgn./srorig;
        dt=1./Fs;
        oldAx=[0:dt:tm]; oldAx(end)=[];
        dtrs1(g,:)=resample(signal(g,:),oldAx,mFs);
    end
    
    % Plot PSDs of the data
    % [pxx1,f] = pwelch(dtrs1(1,:),mFs,0,[1:100],mFs);
    % [pxx2,f] = pwelch(dtrs1(2,:),mFs,0,[1:100],mFs);
    % figure; plot(f,pxx1); hold on; plot(f,pxx2); box off
    % legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
    
    % Find and clean up the PCS data
%     figure; plot(dtrs1(1,:));
    % PCSstr=input('Input timing of the final pulse   ');
      
    PCSstr=SyncDetails(s,1);
    
    dtrst1=dtrs1(:,PCSstr:end);
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,3});
    load(fname);

    mlbltm1=syncR1.strtm;
    % Find and sync the Biopacs / Matlab data.
    figure; plot(syncR1.data(:,3));
    % Biopstr=input('Input timing of the final pulse');
    Biopstr1=SyncDetails(s,2);
    
    trls1=result.data;
    
    % Create Fieldtrip structure
    data=struct;
    % data.fsample=mFs;
    data.label{1,1}='LFP';
    data.label{2,1}='OFC';
    
    ad1=Biopstr1./500;
    % Have checked the timings of this and it works now.
    PCSst1=mlbltm1+ad1;
    
    % Find where a button was pressed %
    trialtype=[];
    for g=1:75;
        if ~isempty(trls1(g).key)
        trialtype(g)=trls1(g).key==37 | trls1(g).key==39;
        end
    end
    whrtrl=find(trialtype==1);
    
    cn=1;
    for g=whrtrl;
        effort(cn,1)=trls1(g).effortIx;
        stake(cn,1)=trls1(g).stakeIx;
        trialstr(cn,1)=trls1(g).startStim-PCSst1;  %1.5 seconds until motor mapping revealed.
        
%         if trls1(g).Yestrial ==1
%             choicestr(g,1)=trls1(g).YesChoice-PCSst1;
%         else
%             choicestr(g,1)=trls1(g).NoChoice-PCSst1;
%         end
%         
        decision(cn,1)=trls1(g).Yestrial;
        st=round(trialstr(cn)*1000)-3000;
        ed=st+6000;
        data.trial{cn}=dtrst1(:,st:ed);
        
%         st=round(choicestr(g)*1000)-3000;
%         ed=st+6000;
        
        % Remember this is dependent on sample rate and should be in seconds.
        data.time{cn}=([st:ed]-st-3000)./1000;
        data.sampleinfo(cn,:)=[st ed];
        %     data.trialinfo(g,1)=trls1(g).effortIx;
        %     data.trialinfo(g,2)=trls1(g).stakeIx;
        cn=cn+1;
    end
    
    %%
    clear syncR1 syncR2 result signal
    close all
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,6});
    load(fname);
       
        % PFC01 - put it back into signal.
    if ~exist('signal','var')
        signal(1,:)=datar(:,1)';        % 1 +2 -0      % 2 = +1 -3
        signal(2,:)=datar(:,3)';        % 3 = +9 -8    % 4 = +11 -10
    end
    
    
    % Synchronise data
    % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
    % Biopacs
    mFs=1000;
    dtrs2=[];
    for g=1:2
        lgn=size(signal,2);
        srorig=Fs;
        tm=lgn./srorig;
        dt=1./Fs;
        oldAx=[0:dt:tm]; oldAx(end)=[];
        dtrs2(g,:)=resample(signal(g,:),oldAx,mFs);
    end
    
    % Plot PSDs of the data
    % [pxx1,f] = pwelch(dtrs1(1,:),mFs,0,[1:100],mFs);
    % [pxx2,f] = pwelch(dtrs1(2,:),mFs,0,[1:100],mFs);
    % figure; plot(f,pxx1); hold on; plot(f,pxx2); box off
    % legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
    
    % Find and clean up the PCS data
    figure; plot(dtrs2(1,:));
    % PCSstr=input('Input timing of the final pulse   ');
    % Note that for this subject - using the first pulse as the end is missing. (Ignore the first
    % deflection which is 300ms before the first pulse, this seems to be a
    % turing on artefact and is not visible in the EMG pulse')
    
    PCSstr=SyncDetails(s,5);
    
    dtrst2=dtrs2(:,PCSstr:end);
    
    % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
    
    fname=fullfile(cd,'DataFiles',FileDetails{s,7});
    load(fname);
    
    mlbltm2=syncR1.strtm;
    % Find and sync the Biopacs / Matlab data.
    figure; plot(syncR1.data(:,3));
    % Biopstr=input('Input timing of the final pulse');
    % Here using the first pulse but be careful to ignore the ECG artefact %
    
    Biopstr2=SyncDetails(s,6);
    
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
    % Find where a button was pressed %
    trialtype=[];
    for g=1:75;
        if ~isempty(trls2(g).key)
        trialtype(g)=trls2(g).key==37 | trls2(g).key==39;
        end
    end
    whrtrl=find(trialtype==1);
    
    cn=1;
    for g=whrtrl;
        effort2(cn,1)=trls2(g).effortIx;
        stake2(cn,1)=trls2(g).stakeIx;
        trialstr2(cn,1)=trls2(g).startStim-PCSst2;  %1.5 seconds until motor mapping revealed.
        
%         if trls2(g).Yestrial ==1
%             choicestr2(g,1)=trls2(g).YesChoice-PCSst2;
%         else
%             choicestr2(g,1)=trls2(g).NoChoice-PCSst2;
%         end
        
        decision2(cn,1)=trls2(g).Yestrial;
        
        st=round(trialstr2(cn)*1000)-3000;
        ed=st+6000;
        data2.trial{cn}=dtrst2(:,st:ed);
        
%         st=round(choicestr2(g)*1000)-3000;
%         ed=st+6000;
%         data2.trial{g}=dtrst2(:,st:ed);
        
        % Remember this is dependent on sample rate and should be in seconds.
        data2.time{cn}=([st:ed]-st-3000)./1000;
        data2.sampleinfo(cn,:)=[st ed];
        %     data2.trialinfo(g,1)=trls2(g).effortIx;
        %     data2.trialinfo(g,2)=trls2(g).stakeIx;
        
        cn=cn+1;
    end
    
    close all;
    clear syncR1 syncR2 result signal
    
    %% Combine together
    
        effortS=effort;
    effortS2=effort2;
    effortS(end)=[];
    effortS=[nan;effortS];
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

    effortS=[effortS; effortS2];
    decisionS=[decisionS;decisionS2];
    stakeS=[stakeS;stakeS2];
    
    % Just add data2 to data 55
    effort=[effort;effort2];
    stake=[stake; stake2];
    decision=[decision;decision2];
    
    blck1length=size(data.trial,2);
    for g=1:size(effort2,1);
        data.trial{g+blck1length}=data2.trial{g};
        data.time{g+blck1length}=data2.time{g};
%         data.trialbl{g+blck1length}=data2.trialbl{75};
    end
    
    % Find where the was a missing trial (not completed).
    whrempty=find(effort==0);
    
    effort(whrempty)=[];
    stake(whrempty)=[];
    decision(whrempty)=[];
    data.trial(whrempty)=[];
    data.time(whrempty)=[];
%     data.trialbl(whrempty)=[];
    
    % data.sampleinfo(76:150,:)=data2.sampleinfo;
    
    for g=1:size(data.trial,2)
        tmp=data.trial{g};
        % Just go off the OFC signal %
        tmp=tmp(2,:);
        trialstd(g,:)=std(tmp,1,2);
        trialmean(g,:)=mean(abs(tmp),2);
    end
    
    mnstd=mean(trialstd);
    trialstdN=trialstd./mnstd;
    
    % figure;
    % subplot(2,1,1); plot(trialstdN(:,1));
    % subplot(2,1,2); plot(trialmean);
    
    whrc=[];
    for g=1:size(trialstdN,2);
        whr=find(trialstdN(:,g)>3);
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
    AllData.desp=strcat(FileDetails{s,1},'Stimulus_Locked');
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