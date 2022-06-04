%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_colors = distinguishable_colors(length(SBJs));

% Analysis parameters:
theta_lim  = [4 7];
beta_lim   = [13 30];    
sbj_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
betapk_bw = 4;
betapk_lim = nan(length(SBJs),2);
for s = 1:length(sbj_beta_pk)
    betapk_lim(s,:) = [sbj_beta_pk(s)-betapk_bw/2 sbj_beta_pk(s)+betapk_bw/2];
end

% Plotting parameters
sem_alpha  = 0.5;
plot_time_lim = [-0.5 2];
stim_time_lim = [0 2];
plot_freq_lim = [2 30];

save_fig = 0;

%% Load data
% clc
% close all
% clear all

% restoredefaultpath
% addpath('C:\BackupRepo\fieldtrip-20190705');
% addpath('C:\Users\slittle\Box\Research\Resources & References\Software\Matlab\General_Code');
% ft_defaults

DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';

[numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);
SyncDetails = numbers;

% Load data files
DataStr=struct;
for s=1:size(FileDetails,1)-1
    
    load(strcat(DataStorage2,'/',FileDetails{s,1},'Stimulus_Locked.mat'));
    DataStr(s).data=AllData;
end

%% Extract group-level power in OFC and subcortical from all SBJ
frqs = DataStr(1).data.TFbl.freq;   % 2:80 Hz
tmes = DataStr(1).data.TFbl.time;   % -3:3 s
ch_lab = DataStr(1).data.TFbl.label;
ofc_ch_ix = find(strcmp(DataStr(1).data.TFbl.label,'OFC'));
lfp_ch_ix = find(strcmp(DataStr(1).data.TFbl.label,'LFP'));

% find frequency indices
for i = 1:2
    [~,theta_freq_ix(i)] = min(abs(frqs-theta_lim(i)));
    [~,beta_freq_ix(i)] = min(abs(frqs-beta_lim(i)));
    for s = 1:length(SBJs)
        [~,betapk_freq_ix(s,i)] = min(abs(frqs-betapk_lim(s,i)));
    end
end

% Compute mean theta and beta power across trials in PFC
theta_ofc   = nan(length(SBJs),length(tmes));
beta_ofc    = nan(length(SBJs),length(tmes));
betapk_ofc  = nan(length(SBJs),length(tmes));
pwr_ofc     = nan(length(SBJs),length(frqs),length(tmes));
theta_lfp   = nan(length(SBJs),length(tmes));
beta_lfp    = nan(length(SBJs),length(tmes));
betapk_lfp  = nan(length(SBJs),length(tmes));
pwr_lfp     = nan(length(SBJs),length(frqs),length(tmes));
for s=1:length(SBJs)
    % Check data matches across SBJs
    if ~all(frqs==DataStr(s).data.TFbl.freq); error('frequencies mismatch'); end
    if ~all(tmes==DataStr(s).data.TFbl.time); error('times mismatch'); end
    if ~strcmp(DataStr(s).data.TFbl.label{ofc_ch_ix},'OFC'); error('OFC channel index mismatch'); end
    if ~strcmp(DataStr(s).data.TFbl.label{lfp_ch_ix},'LFP'); error('LFP channel index mismatch'); end
    
    % avg across trials using baseline corrected .TFbl, not raw .TF
    pwr_mean = squeeze(mean(DataStr(s).data.TFbl.powspctrm,1));
    
    theta_ofc(s,:)  = squeeze(mean(pwr_mean(ofc_ch_ix,theta_freq_ix(1):theta_freq_ix(2),:),2));
    beta_ofc(s,:)   = squeeze(mean(pwr_mean(ofc_ch_ix,beta_freq_ix(1):beta_freq_ix(2),:),2));
    betapk_ofc(s,:) = squeeze(mean(pwr_mean(ofc_ch_ix,betapk_freq_ix(s,1):betapk_freq_ix(s,2),:),2));
    pwr_ofc(s,:,:)  = squeeze(pwr_mean(ofc_ch_ix,:,:)); % [sbj, freq, time]
    
    theta_lfp(s,:)  = squeeze(mean(pwr_mean(lfp_ch_ix,theta_freq_ix(1):theta_freq_ix(2),:),2));
    beta_lfp(s,:)   = squeeze(mean(pwr_mean(lfp_ch_ix,beta_freq_ix(1):beta_freq_ix(2),:),2));
    betapk_lfp(s,:) = squeeze(mean(pwr_mean(lfp_ch_ix,betapk_freq_ix(s,1):betapk_freq_ix(s,2),:),2));
    pwr_lfp(s,:,:) = squeeze(pwr_mean(lfp_ch_ix,:,:)); % [sbj, freq, time]
end

%% Plot Group-level OFC and LFP power
% find time indices
[~,stim_start_ix] = min(abs(tmes-stim_time_lim(1)));
[~,stim_end_ix]   = min(abs(tmes-stim_time_lim(2)));
[~,tme_start_ix] = min(abs(tmes-plot_time_lim(1)));
[~,tme_end_ix]   = min(abs(tmes-plot_time_lim(2)));
plt_tme_idx = tme_start_ix:tme_end_ix;
plot_tmes = tmes(plt_tme_idx);
[~,stim_ix] = min(abs(plot_tmes-0));

% Compute standard error
theta_ofc_sem = std(theta_ofc,[],1)./sqrt(length(SBJs));
theta_lfp_sem = std(theta_lfp,[],1)./sqrt(length(SBJs));
beta_ofc_sem  = std(beta_ofc,[],1)./sqrt(length(SBJs));
beta_lfp_sem  = std(beta_lfp,[],1)./sqrt(length(SBJs));

figure('units','norm','OuterPosition',[0 0 1 1]);

% ================================ PLOT OFC ================================
% Plot OFC theta
subplot(2,3,1);
shadedErrorBar(plot_tmes, mean(theta_ofc(:,plt_tme_idx),1),theta_ofc_sem(:,plt_tme_idx));%, 'transparent',sem_alpha);
line([0 0],ylim,'Color','k','LineStyle','-');
xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
ylabel('OFC Theta Power');
title(['Group averaged (n=' num2str(length(SBJs)) ') OFC theta power']);
% xlim([-0.5 2]);
set(gca,'FontSize',16);

% Plot OFC beta
subplot(2,3,2); hold on;
b_line   = shadedErrorBar(plot_tmes, mean(beta_ofc(:,plt_tme_idx),1),...
    beta_ofc_sem(:,plt_tme_idx),'lineprops',{'Color','r'});%, 'transparent',sem_alpha);
bpk_line = shadedErrorBar(plot_tmes, mean(betapk_ofc(:,plt_tme_idx),1),...
    beta_ofc_sem(:,plt_tme_idx),'lineprops',{'Color','b'});%, 'transparent',sem_alpha);
line([0 0],ylim,'Color','k','LineStyle','-');

xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
ylabel('OFC Beta Power');
title(['Group averaged (n=' num2str(length(SBJs)) ') OFC beta power']);
% xlim([-0.5 2]);
set(gca,'FontSize',16);
legend([b_line.mainLine, bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ' Hz)'],...
    'SBJ-specific beta'},'Location','best');

% Plot OFC TFR
subplot(2,3,3);
% clims=[0 3.5];
imagesc(plot_tmes, frqs, squeeze(mean(pwr_ofc(:,:,plt_tme_idx),1)));%,clims);
line([0 0],ylim,'Color','k','LineStyle','-');

set(gca,'YDir','normal')
% xlim([0 1.5]);
xlabel('Time (s)');
ylim(plot_freq_lim);
ylabel('Frequency (Hz)');
colorbar;
title(['Group averaged (n=' num2str(length(SBJs)) ') OFC TFR']);
set(gca,'FontSize',16);

% ================================ PLOT LFP ================================
% Plot LFP theta
subplot(2,3,4);
shadedErrorBar(plot_tmes, mean(theta_lfp(:,plt_tme_idx),1),...
    theta_lfp_sem(:,plt_tme_idx));%, 'transparent',sem_alpha);
line([0 0],ylim,'Color','k','LineStyle','-');
xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
ylabel('LFP Theta Power');
title(['Group averaged (n=' num2str(length(SBJs)) ') LFP theta power']);
% xlim([-0.5 2]);
set(gca,'FontSize',16);

% Plot LFP beta
subplot(2,3,5); hold on;
b_line   = shadedErrorBar(plot_tmes, mean(beta_lfp(:,plt_tme_idx),1),...
    beta_lfp_sem(:,plt_tme_idx),'lineprops',{'Color','r'});%, 'transparent',sem_alpha);
bpk_line = shadedErrorBar(plot_tmes, mean(betapk_ofc(:,plt_tme_idx),1),...
    beta_ofc_sem(:,plt_tme_idx),'lineprops',{'Color','b'});%, 'transparent',sem_alpha);
line([0 0],ylim,'Color','k','LineStyle','-');

xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
ylabel('LFP Beta Power');
title(['Group averaged (n=' num2str(length(SBJs)) ') LFP beta power']);
% xlim([-0.5 2]);
set(gca,'FontSize',16);
legend([b_line.mainLine, bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ' Hz)'],...
    'SBJ-specific beta'},'Location','best');

% Plot LFP TFR
subplot(2,3,6);
% clims=[0 3.5];
imagesc(plot_tmes, frqs, squeeze(mean(pwr_lfp(:,:,plt_tme_idx),1)));%,clims);
line([0 0],ylim,'Color','k','LineStyle','-');

set(gca,'YDir','normal')
% xlim([0 1.5]);
xlabel('Time (s)');
ylim(plot_freq_lim);
ylabel('Frequency (Hz)');
colorbar;
title(['Group averaged (n=' num2str(length(SBJs)) ') LFP TFR']);
set(gca,'FontSize',16);

%% Extract SBJ-level power in OFC and LFP
% Compute SEM across trials
sbj_tfr         = nan(length(SBJs),length(ch_lab),length(frqs),length(tmes));
sbj_beta_max    = nan(length(SBJs),length(ch_lab));
sbj_betad_max   = nan(length(SBJs),length(ch_lab));
sbj_theta_sem   = nan(length(SBJs),length(ch_lab),length(tmes));
sbj_beta_sem    = nan(length(SBJs),length(ch_lab),length(tmes));
sbj_betapk_sem  = nan(length(SBJs),length(ch_lab),length(tmes));
for s=1:length(SBJs)
    sbj_pwr = DataStr(s).data.TFbl.powspctrm;
    sbj_tfr(s,:,:,:) = squeeze(nanmean(sbj_pwr,1));
    psd = squeeze(nanmean(sbj_tfr(s,:,:,:),4));
    frqs_beta = frqs(beta_freq_ix(1):beta_freq_ix(2));
    [~,max_ix] = max(psd(:,beta_freq_ix(1):beta_freq_ix(2)),[],2);
    [~,maxd_ix] = max(diff(psd(:,beta_freq_ix(1):beta_freq_ix(2))),[],2);
    
    sbj_beta_max(s,:) = frqs_beta(max_ix);
    sbj_betad_max(s,:) = frqs_beta(maxd_ix);
    
    % Average within frequency band
    sbj_theta  = squeeze(mean(sbj_pwr(:,:,theta_freq_ix(1):theta_freq_ix(2),:),3));
    sbj_beta   = squeeze(mean(sbj_pwr(:,:,beta_freq_ix(1):beta_freq_ix(2),:),3));
    sbj_betapk = squeeze(mean(sbj_pwr(:,:,betapk_freq_ix(s,1):betapk_freq_ix(s,2),:),3));
    
    sbj_theta_sem(s,:,:)  = squeeze(std(sbj_theta,[],1))./sqrt(size(sbj_theta,1));
    sbj_beta_sem(s,:,:)   = squeeze(std(sbj_beta,[],1))./sqrt(size(sbj_beta,1));
    sbj_betapk_sem(s,:,:) = squeeze(std(sbj_betapk,[],1))./sqrt(size(sbj_betapk,1));
    
%     sbj_theta_lfp_sem(s,:)  = squeeze(std(sbj_pwr(:,lfp_ch_ix,theta_freq_ix(1):theta_freq_ix(2),:),[],1))./sqrt(size(sbj_pwr,1));
%     sbj_beta_lfp_sem(s,:)   = squeeze(std(sbj_pwr(:,lfp_ch_ix,beta_freq_ix(1):beta_freq_ix(2),:),[],1))./sqrt(size(sbj_pwr,1));
%     sbj_betapk_lfp_sem(s,:) = squeeze(std(sbj_pwr(:,lfp_ch_ix,betapk_freq_ix(1):betapk_freq_ix(2),:),[],1))./sqrt(size(sbj_pwr,1));
end


%% Plot SBJ-level OFC and LPF power
for ch_ix = 1:length(ch_lab)
    fig_name = ['SBJ power for ' ch_lab{ch_ix}];
    figure('Name',fig_name,'units','norm','OuterPosition',[0 0 1 1]);
    for s = 1:length(SBJs)
        % Plot TFR
        subplot(length(SBJs),4,(s*4)-3);% hold on;
        imagesc(plot_tmes, frqs, squeeze(sbj_tfr(s,ch_ix,:,plt_tme_idx)));%,clims);
        line([0 0],ylim,'Color','k','LineStyle','-');
        set(gca,'YDir','normal')
        % xlim([0 1.5]);
        xlabel('Time (s)');
        ylim(plot_freq_lim);
        ylabel('Frequency (Hz)');
        colorbar;
        title([SBJs{s} ': ' ch_lab{ch_ix} ' TFR']);
        set(gca,'FontSize',16);
        
        % Plot PSD
        subplot(length(SBJs),4,(s*4)-2); hold on;
        stim_psd_line = plot(frqs,squeeze(nanmean(sbj_tfr(s,ch_ix,:,stim_start_ix:stim_end_ix),4)),'r');
        psd_line = plot(frqs,squeeze(nanmean(sbj_tfr(s,ch_ix,:,:),4)),'b');
%         line([sbj_beta_max(s,ch_ix) sbj_beta_max(s,ch_ix)],ylim,'Color','b');
        psdd_line = plot(frqs(2:end),diff(squeeze(nanmean(sbj_tfr(s,ch_ix,:,:),4))),'k');
%         line([sbj_betad_max(s,ch_ix) sbj_betad_max(s,ch_ix)],ylim,'Color','k');
        xlabel('Frequency (Hz)');
        ylabel('Power (norm)');%'Power mV^2/Hz');
        xlim(plot_freq_lim);
        legend([psd_line stim_psd_line psdd_line],{'PSD',...
            ['PSD (' num2str(stim_time_lim(1)),'-',num2str(stim_time_lim(2)),' s)'],'diff(PSD)'},'Location','best');
        title([SBJs{s} ': ' ch_lab{ch_ix} ' PSD']);
        set(gca,'FontSize',16);
        
        % Plot theta
        subplot(length(SBJs),4,(s*4)-1); hold on;
        shadedErrorBar(plot_tmes, theta_lfp(s,plt_tme_idx),sbj_theta_sem(s,ch_ix,plt_tme_idx));
        line([0 0],ylim,'Color','k','LineStyle','-');
        xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
        ylabel('Theta Power');
        title([SBJs{s} ': ' ch_lab{ch_ix} ' Theta']);
        % xlim([-0.5 2]);
        set(gca,'FontSize',16);
        
        % Plot beta
        subplot(length(SBJs),4,s*4); hold on;
        shadedErrorBar(plot_tmes, beta_lfp(s,plt_tme_idx),sbj_beta_sem(s,ch_ix,plt_tme_idx),'lineprops',{'Color','r'});
        shadedErrorBar(plot_tmes, betapk_lfp(s,plt_tme_idx),sbj_betapk_sem(s,ch_ix,plt_tme_idx),'lineprops',{'Color','b'});
        line([0 0],ylim,'Color','k','LineStyle','-');
        xlabel('Time (s)','Fontsize',10,'Fontname','Arial');
        ylabel('Beta Power');
        title([SBJs{s} ': ' ch_lab{ch_ix} ' Beta']);
        % xlim([-0.5 2]);
        set(gca,'FontSize',16);
        legend([num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ' Hz'],...
            [num2str(betapk_lim(s,1)) '-' num2str(betapk_lim(s,2)) ' Hz'],'Location','best');
    end
end

%% Plot PFC01 singel-trial stacks for outliers ~1.5s
s = find(strcmp(SBJs,'PFC01'));
% Get PFC01 single-trial theta/beta
sbj_pwr = DataStr(s).data.TFbl.powspctrm;
sbj_theta  = squeeze(mean(sbj_pwr(:,:,theta_freq_ix(1):theta_freq_ix(2),:),3));
sbj_beta   = squeeze(mean(sbj_pwr(:,:,beta_freq_ix(1):beta_freq_ix(2),:),3));

figure('Name',['PFC01 single-trial stacks']);
for ch_ix = 1:length(ch_lab)
    % Plot Theta single-trial stack
    subplot(2,2,ch_ix*2-1);
    imagesc(plot_tmes, 1:size(sbj_theta,1), squeeze(sbj_theta(:,ch_ix,plt_tme_idx)));%,clims);
    line([0 0],ylim,'Color','k','LineStyle','-');
    set(gca,'YDir','normal')
    % xlim([0 1.5]);
    xlabel('Time (s)');
    ylabel('Trials');
    colorbar;
    title([SBJs{s} ': ' ch_lab{ch_ix} ' single-trial theta']);
    set(gca,'FontSize',16);
    
    % Plot Beta single-trial stack
    subplot(2,2,ch_ix*2);
    imagesc(plot_tmes, 1:size(sbj_beta,1), squeeze(sbj_beta(:,ch_ix,plt_tme_idx)));%,clims);
    line([0 0],ylim,'Color','k','LineStyle','-');
    set(gca,'YDir','normal')
    % xlim([0 1.5]);
    xlabel('Time (s)');
    ylabel('Trials');
    colorbar;
    title([SBJs{s} ': ' ch_lab{ch_ix} ' single-trial beta']);
    set(gca,'FontSize',16);
end

%% Plot problematic trials
outlier_sd_thresh = 3;
sbj_theta_pk_trl = max(sbj_theta,[],3);
mn_pk = mean(sbj_theta_pk_trl,1);
sd_pk = std(sbj_theta_pk_trl,[],1);
bad_trl_ix = [];
for ch_ix = 1:2
    bad_trl_ix = [bad_trl_ix; find(sbj_theta_pk_trl(:,ch_ix)>mn_pk(ch_ix)+(sd_pk(ch_ix)*outlier_sd_thresh))];
end
for t = 1:length(bad_trl_ix)
    fprintf('%d = %.01f OFC and %.01f LFP\n',bad_trl_ix(t),sbj_theta_pk_trl(bad_trl_ix(t),1),sbj_theta_pk_trl(bad_trl_ix(t),2));
end

cfg_reject = [];
cfg_reject.method = 'summary';
ft_rejectvisual(cfg_reject,DataStr(s).data.raw);

cfg_plot = [];
cfg_plot.viewmode = 'vertical';
plt_data = DataStr(s).data.raw;
plt_data = rmfield(plt_data,'sampleinfo');
ft_databrowser(cfg_plot,plt_data);

