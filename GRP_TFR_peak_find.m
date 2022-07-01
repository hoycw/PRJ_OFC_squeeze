%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc

addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_id = 'TFRw_S25t2_noBsln_fl2t40_c7';

if contains(an_id,'_S')
    psd_win_lim = [0 2];
    an_lim = [0.5 1];
elseif contains(an_id,'_D')
    psd_win_lim = [-0.5 0];
    an_lim = [-0.5 0];
end

% Analysis parameters:
theta_lim  = [3 9];
beta_lim   = [13 30];    
sbj_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
thetapk_bw = 4;
betapk_bw = 4;

% Plotting parameters
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

freq_ticks = 5:5:35;
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
fig_dir   = [prj_dir 'results/TFR/' an_id '/pk_find/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load TFR data
theta_pk = nan([length(SBJs) 2]);
beta_pk = nan([length(SBJs) 2]);
thetapk_lim = nan(length(SBJs),2,2);
betapk_lim  = nan(length(SBJs),2,2);
for s = 1:length(SBJs)
    %% Load data
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    
    % Initialize data
    if s==1
        theta_avg  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        thetapk_avg= nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_avg   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        betapk_avg = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        theta_sem  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        thetapk_sem= nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_sem   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        betapk_sem = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        time_vec = tmp.tfr.time;
        freq_vec = tmp.tfr.freq;
        tfr      = nan([length(SBJs) 2 numel(tmp.tfr.freq) numel(tmp.tfr.time)]);
        psd      = nan([length(SBJs) 2 numel(tmp.tfr.freq)]);
    end
    
    % Average across trials
    tfr(s,:,:,:) = squeeze(nanmean(tmp.tfr.powspctrm,1));
    
    % Compute PSD
    cfgs = [];
    cfgs.latency = psd_win_lim;
    psd_tfr = ft_selectdata(cfgs,tmp.tfr);
    psd(s,:,:) = squeeze(nanmean(nanmean(psd_tfr.powspctrm,1),4));
    
    % Find PSD peaks
    for ch_ix = 1:2
        % Theta
        [~,min_ix] = min(abs(freq_vec-theta_lim(1)));
        [~,max_ix] = min(abs(freq_vec-theta_lim(2)));
        [pk,pk_ix] = findpeaks(squeeze(psd(s,ch_ix,min_ix:max_ix)));
        if numel(pk_ix)~=0
            if numel(pk_ix)>1
                [~,best_ix] = max(pk);
                pk_ix = pk_ix(best_ix);
                warning('Warning: %s has more than one peak in %s theta! taking the max...',SBJs{s},tmp.tfr.label{ch_ix});
            end
            theta_pk(s,ch_ix) = freq_vec(min_ix+pk_ix-1);
        else
            warning('Warning: No peaks found for %s in theta band...');
        end
        
        % Beta
        [~,min_ix] = min(abs(freq_vec-beta_lim(1)));
        [~,max_ix] = min(abs(freq_vec-beta_lim(2)));
        [pk,pk_ix] = findpeaks(squeeze(psd(s,ch_ix,min_ix:max_ix)));
        if numel(pk_ix)~=0
            if numel(pk_ix)>1
                [~,best_ix] = max(pk);
                pk_ix = pk_ix(best_ix);
                warning('Warning: %s has more than one peak in %s beta! taking the max...',SBJs{s},tmp.tfr.label{ch_ix});
            end
            beta_pk(s,ch_ix) = freq_vec(min_ix+pk_ix-1);
        else
            warning('Warning: No peaks found for %s in theta band...');
        end
        
        % Compute bandwidth
        if ~isnan(theta_pk(s,ch_ix))
            thetapk_lim(s,ch_ix,:) = [theta_pk(s,ch_ix,1)-thetapk_bw/2 theta_pk(s,ch_ix,1)+thetapk_bw/2];
        else
            thetapk_lim(s,ch_ix,:) = theta_lim;
        end
        if ~isnan(beta_pk(s,ch_ix))
            betapk_lim(s,ch_ix,:) = [beta_pk(s,ch_ix,1)-betapk_bw/2 beta_pk(s,ch_ix,1)+betapk_bw/2];
        else
            betapk_lim(s,ch_ix,:) = beta_lim;
        end
        
        % Compute power bands
        cfg = [];
        cfg.channel = tmp.tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.frequency = theta_lim;
        theta_pow = ft_selectdata(cfg, tmp.tfr);
        theta_sem(s,ch_ix,:) = squeeze(nanstd(theta_pow.powspctrm,[],1))./sqrt(size(theta_pow.powspctrm,1));
        theta_avg(s,ch_ix,:) = squeeze(nanmean(theta_pow.powspctrm,1));
        
        cfg.frequency = squeeze(thetapk_lim(s,ch_ix,:))';
        thetapk_pow = ft_selectdata(cfg, tmp.tfr);
        thetapk_sem(s,ch_ix,:) = squeeze(nanstd(thetapk_pow.powspctrm,[],1))./sqrt(size(thetapk_pow.powspctrm,1));
        thetapk_avg(s,ch_ix,:) = squeeze(nanmean(thetapk_pow.powspctrm,1));
        
        cfg.frequency = beta_lim;
        beta_pow = ft_selectdata(cfg, tmp.tfr);
        beta_sem(s,ch_ix,:) = squeeze(nanstd(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
        beta_avg(s,ch_ix,:) = squeeze(nanmean(beta_pow.powspctrm,1));
        
        cfg.frequency = squeeze(betapk_lim(s,ch_ix,:))';
        betapk_pow = ft_selectdata(cfg, tmp.tfr);
        betapk_sem(s,ch_ix,:) = squeeze(nanstd(betapk_pow.powspctrm,[],1))./sqrt(size(betapk_pow.powspctrm,1));
        betapk_avg(s,ch_ix,:) = squeeze(nanmean(betapk_pow.powspctrm,1));
    end
end

%% Plot SBJ TFRs
% Get frequency and time ticks
freq_tick_ix = nan(size(freq_ticks));
for f = 1:numel(freq_ticks)
    [~,freq_tick_ix(f)] = min(abs(freq_vec-freq_ticks(f)));
end
time_ticks = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
time_tick_ix = nan(size(time_ticks));
for t = 1:numel(time_ticks)
    [~,time_tick_ix(t)] = min(abs(time_vec-time_ticks(t)));
end
if ~any(time_ticks==0); error('cant plot event with no tick at 0'); end

for s = 1:length(SBJs)
    %% Plot TFRs for each SBJ
    fig_name = [SBJs{s} '_TFR_theta_beta_' an_id];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        subplot(2,4,ch_ix*4-3); hold on;
        % Get color lims per condition
        vals = tfr(s,ch_ix,:,:);
        clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
        
        % Plot TFR
        imagesc(squeeze(tfr(s,ch_ix,:,:)));
        caxis(clim);
        colorbar;
        
        % Plot Events
        line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        if ch_ix==1; ch_lab = 'LFP'; else ch_lab = sbj_pfc_roi{s}; end
        title([ch_lab '- ' an_id], 'interpreter', 'none');
        set(gca,'XLim',[1 numel(time_vec)]);
        set(gca,'XTick', time_tick_ix);
        set(gca,'XTickLabels', time_ticks);
        set(gca,'YLim',[1 numel(freq_vec)]);
        set(gca,'YTick',freq_tick_ix);
        set(gca,'YTickLabels',freq_ticks);
        set(gca,'YDir','normal');
%         ylim([min(freq_vec) max(freq_vec)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        set(gca,'FontSize',16);
        
        %% Plot PSD
        subplot(2,4,ch_ix*4-2); hold on;
        lfp_line = plot(freq_vec,squeeze(psd(s,1,:)),'Color','b');
        pfc_line = plot(freq_vec,squeeze(psd(s,2,:)),'Color',[255,165,0]./255); %orange
        ylims = ylim;
        if ~isnan(theta_pk(s,ch_ix))
            patch([thetapk_lim(s,ch_ix,1) thetapk_lim(s,ch_ix,1) ...
                thetapk_lim(s,ch_ix,2) thetapk_lim(s,ch_ix,2)],...
                [ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0.1);
        end
        if ~isnan(beta_pk(s,ch_ix))
            patch([betapk_lim(s,ch_ix,1) betapk_lim(s,ch_ix,1) ...
                betapk_lim(s,ch_ix,2) betapk_lim(s,ch_ix,2)],...
                [ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0.1);
        end
        xlabel('Frequency (Hz)');
        ylabel('Power (norm)');
        legend([lfp_line pfc_line],{'LFP',sbj_pfc_roi{s}},'Location','best');
        title(['PSD (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
        set(gca,'FontSize',16);
        
        %% Plot theta
        subplot(2,4,ch_ix*4-1); hold on;
        
        t_line = shadedErrorBar(time_vec, squeeze(theta_avg(s,ch_ix,:)),...
            squeeze(theta_sem(s,ch_ix,:)),'lineprops',{'Color','r'});
        tpk_line = shadedErrorBar(time_vec, squeeze(thetapk_avg(s,ch_ix,:)),...
            squeeze(thetapk_sem(s,ch_ix,:)),'lineprops',{'Color','b'});
        line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_lab ' theta'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend([t_line.mainLine tpk_line.mainLine],{['theta (' num2str(theta_lim(1)) '-' num2str(theta_lim(2)) ...
            ' Hz)'],['SBJ theta (' num2str(thetapk_lim(s,ch_ix,1)) '-' num2str(thetapk_lim(s,ch_ix,2)) ' Hz)']},'Location','best');
        set(gca,'FontSize',16);
        
        %% Plot beta
        subplot(2,4,ch_ix*4); hold on;
        
        b_line = shadedErrorBar(time_vec, squeeze(beta_avg(s,ch_ix,:)),...
            squeeze(beta_sem(s,ch_ix,:)),'lineprops',{'Color','r'});
        bpk_line = shadedErrorBar(time_vec, squeeze(betapk_avg(s,ch_ix,:)),...
            squeeze(betapk_sem(s,ch_ix,:)),'lineprops',{'Color','b'});
        line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_lab ' beta'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
            ' Hz)'],['SBJ beta (' num2str(betapk_lim(s,ch_ix,1)) '-' num2str(betapk_lim(s,ch_ix,2)) ' Hz)']},'Location','best');
        set(gca,'FontSize',16);
        
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Average at group level
lfp_tfr_grp        = squeeze(nanmean(tfr(:,1,:,:),1));
lfp_theta_grp_avg  = squeeze(mean(theta_avg(:,1,:),1));
lfp_theta_grp_sem  = squeeze(std(theta_avg(:,1,:),[],1)./sqrt(size(theta_avg(:,1,:),1)));
lfp_thetapk_grp_avg  = squeeze(mean(thetapk_avg(:,1,:),1));
lfp_thetapk_grp_sem  = squeeze(std(thetapk_avg(:,1,:),[],1)./sqrt(size(thetapk_avg(:,1,:),1)));
lfp_beta_grp_avg   = squeeze(mean(beta_avg(:,1,:),1));
lfp_beta_grp_sem   = squeeze(std(beta_avg(:,1,:),[],1)./sqrt(size(beta_avg(:,1,:),1)));
lfp_betapk_grp_avg = squeeze(mean(betapk_avg(:,1,:),1));
lfp_betapk_grp_sem = squeeze(std(betapk_avg(:,1,:),[],1)./sqrt(size(betapk_avg(:,1,:),1)));

fpc_tfr_grp        = squeeze(nanmean(tfr(strcmp(sbj_pfc_roi,'FPC'),2,:,:),1));
fpc_theta_grp_avg  = squeeze(mean(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_theta_grp_sem  = squeeze(std(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));
fpc_thetapk_grp_avg  = squeeze(mean(thetapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_thetapk_grp_sem  = squeeze(std(thetapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(thetapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));
fpc_beta_grp_avg   = squeeze(mean(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_beta_grp_sem   = squeeze(std(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));
fpc_betapk_grp_avg = squeeze(mean(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_betapk_grp_sem = squeeze(std(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));

ofc_tfr_grp        = squeeze(nanmean(tfr(strcmp(sbj_pfc_roi,'OFC'),2,:,:),1));
ofc_theta_grp_avg  = squeeze(mean(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_theta_grp_sem  = squeeze(std(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));
ofc_thetapk_grp_avg  = squeeze(mean(thetapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_thetapk_grp_sem  = squeeze(std(thetapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(thetapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));
ofc_beta_grp_avg   = squeeze(mean(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_beta_grp_sem   = squeeze(std(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));
ofc_betapk_grp_avg = squeeze(mean(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_betapk_grp_sem = squeeze(std(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));

%% Plot TFRs for the GROUP
fig_name = ['GRP_TFR_theta_beta_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);

%--------------------------------------------------------------------------
%---------------------- LFP ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,4,1); hold on;
% Get color lims per condition
clim = [prctile(lfp_tfr_grp(:),plt.clim_perc(1)) prctile(lfp_tfr_grp(:),plt.clim_perc(2))];

imagesc(lfp_tfr_grp);
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
    'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title(['LFP- ' an_id], 'interpreter', 'none');
set(gca,'XLim',[1 numel(time_vec)]);
set(gca,'XTick', time_tick_ix);
set(gca,'XTickLabels', time_ticks);
set(gca,'YLim',[1 numel(freq_vec)]);
set(gca,'YTick',freq_tick_ix);
set(gca,'YTickLabels',freq_ticks);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot PSD
subplot(3,4,2); hold on;
plot(freq_vec,squeeze(nanmean(psd(:,1,:),1)),'Color','b');
xlabel('Frequency (Hz)');
ylabel('Power (norm)');
title(['GRP LFP PSD (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
set(gca,'FontSize',16);

% Plot theta
subplot(3,4,3); hold on;

t_line = shadedErrorBar(time_vec, lfp_theta_grp_avg, lfp_theta_grp_sem,'lineprops',{'Color','r'});
tpk_line = shadedErrorBar(time_vec, lfp_thetapk_grp_avg, lfp_thetapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('LFP theta', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([t_line.mainLine tpk_line.mainLine],{['theta (' num2str(theta_lim(1)) '-' num2str(theta_lim(2)) ...
    ' Hz)'],['SBJ thWeta (' num2str(mean(thetapk_lim(:,1,1))) '-' num2str(mean(thetapk_lim(:,1,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

% Plot beta
subplot(3,4,4); hold on;

b_line = shadedErrorBar(time_vec, lfp_beta_grp_avg, lfp_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, lfp_betapk_grp_avg, lfp_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('LFP beta', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(mean(betapk_lim(:,1,1))) '-' num2str(mean(betapk_lim(:,1,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

%--------------------------------------------------------------------------
%---------------------- FPC ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,4,5); hold on;
% Get color lims per condition
clim = [prctile(fpc_tfr_grp(:),plt.clim_perc(1)) prctile(fpc_tfr_grp(:),plt.clim_perc(2))];

imagesc(fpc_tfr_grp);
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
    'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title(['FPC- ' an_id], 'interpreter', 'none');
set(gca,'XLim',[1 numel(time_vec)]);
set(gca,'XTick', time_tick_ix);
set(gca,'XTickLabels', time_ticks);
set(gca,'YLim',[1 numel(freq_vec)]);
set(gca,'YTick',freq_tick_ix);
set(gca,'YTickLabels',freq_ticks);xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot PSD
subplot(3,4,6); hold on;
plot(freq_vec,squeeze(nanmean(psd(strcmp(sbj_pfc_roi,'FPC'),2,:),1)),'Color',[255,165,0]./255);
xlabel('Frequency (Hz)');
ylabel('Power (norm)');
title(['GRP FPC PSD (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
set(gca,'FontSize',16);

% Plot theta
subplot(3,4,7); hold on;

t_line = shadedErrorBar(time_vec, fpc_theta_grp_avg, fpc_theta_grp_sem,'lineprops',{'Color','r'});
tpk_line = shadedErrorBar(time_vec, fpc_thetapk_grp_avg, fpc_thetapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('FPC theta', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([t_line.mainLine tpk_line.mainLine],{['theta (' num2str(theta_lim(1)) '-' num2str(theta_lim(2)) ...
    ' Hz)'],['SBJ theta (' num2str(mean(thetapk_lim(strcmp(sbj_pfc_roi,'FPC'),2,1)))...
    '-' num2str(mean(thetapk_lim(strcmp(sbj_pfc_roi,'FPC'),2,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

% Plot beta
subplot(3,4,8); hold on;

b_line = shadedErrorBar(time_vec, fpc_beta_grp_avg, fpc_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, fpc_betapk_grp_avg, fpc_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('FPC beta:', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(mean(betapk_lim(strcmp(sbj_pfc_roi,'FPC'),2,1)))...
    '-' num2str(mean(betapk_lim(strcmp(sbj_pfc_roi,'FPC'),2,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

%--------------------------------------------------------------------------
%---------------------- OFC ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,4,9); hold on;
% Get color lims per condition
clim = [prctile(ofc_tfr_grp(:),plt.clim_perc(1)) prctile(ofc_tfr_grp(:),plt.clim_perc(2))];

imagesc(ofc_tfr_grp);
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
    'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title(['OFC- ' an_id], 'interpreter', 'none');
set(gca,'XLim',[1 numel(time_vec)]);
set(gca,'XTick', time_tick_ix);
set(gca,'XTickLabels', time_ticks);
set(gca,'YLim',[1 numel(freq_vec)]);
set(gca,'YTick',freq_tick_ix);
set(gca,'YTickLabels',freq_ticks);xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot PSD
subplot(3,4,10); hold on;
plot(freq_vec,squeeze(nanmean(psd(strcmp(sbj_pfc_roi,'OFC'),2,:),1)),'Color',[255,165,0]./255);
xlabel('Frequency (Hz)');
ylabel('Power (norm)');
title(['GRP OFC PSD (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
set(gca,'FontSize',16);

% Plot theta
subplot(3,4,11); hold on;

t_line = shadedErrorBar(time_vec, ofc_theta_grp_avg, ofc_theta_grp_sem,'lineprops',{'Color','r'});
tpk_line = shadedErrorBar(time_vec, ofc_thetapk_grp_avg, ofc_thetapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('OFC theta', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([t_line.mainLine tpk_line.mainLine],{['theta (' num2str(theta_lim(1)) '-' num2str(theta_lim(2)) ...
    ' Hz)'],['SBJ theta (' num2str(mean(thetapk_lim(strcmp(sbj_pfc_roi,'OFC'),2,1)))...
    '-' num2str(mean(thetapk_lim(strcmp(sbj_pfc_roi,'OFC'),2,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

% Plot beta
subplot(3,4,12); hold on;

b_line = shadedErrorBar(time_vec, ofc_beta_grp_avg, ofc_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, ofc_betapk_grp_avg, ofc_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title(['OFC beta'], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(mean(betapk_lim(strcmp(sbj_pfc_roi,'OFC'),2,1)))...
    '-' num2str(mean(betapk_lim(strcmp(sbj_pfc_roi,'OFC'),2,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot theta power
% fig_name = ['GRP_theta_' an_id];
% figure('Name',fig_name,'units','normalized',...
%     'outerposition',[0 0 1 1],'Visible',fig_vis);
% 
% for ch_ix = 1:numel(tfr_grp{1}.label)
%     for an_ix = 1:numel(an_ids)
%         subplot(length(tfr_grp{an_ix}.label),length(an_ids),ch_ix*length(an_ids)-an_ix+1);
%         % Get color lims per condition
%         vals = tfr_grp{an_ix}.powspctrm(ch_ix,:,:);
%         clim = [min(vals(:)) max(vals(:))];
%         
%         % Plot TFR
%         imagesc(tfr_grp{an_ix}.time, tfr_grp{an_ix}.freq, squeeze(tfr_grp{an_ix}.powspctrm(ch_ix,:,:)));
%         set(gca,'YDir','normal');
%         caxis(clim);
%         colorbar;
%         
%         % Plot Events
%         line([0 0],ylim,...
%             'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
%             'LineStyle',plt.evnt_styles{1});
%         
%         % Axes and parameters
%         title([tfr_grp{an_ix}.label{ch_ix} '- ' an_ids{an_ix}], 'interpreter', 'none');
%         set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
%         set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
%         xlabel('Time (s)');
%         ylabel('Frequency (Hz)');
%         set(gca,'FontSize',16);
%     end
% end
% 
% % Save Figure
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end
% 

%% Time Frequency analysis on the dataset.
%     cfg=[];
%     cfg.trials          = 'all';
%     cfg.keeptrials      = 'yes';
%     cfg.output          = 'pow';
%     cfg.method          = 'mtmconvol';
%     cfg.taper           = 'hanning';
%     cfg.foi             = 2:80;
%     % data.sampleinfo(1,2)
%     cfg.toi             = -3:0.02:3;
%     cfg.pad             = 'maxperlen';
%     % cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;
%     % Frequecy dependent wavelet length %
%     cfg.t_ftimwin       = 4./cfg.foi;
%     freqout             = ft_freqanalysis(cfg, data);
%
%     cfg.baseline     = [-1 0];
%     cfg.baselinetype = 'zscore';    % 'absolute', 'relative' ratio, 'relchange' %, 'normchange', 'db', 'zscore'.
%     cfg.parameter    = 'powspctrm';
%     freqoutbl = ft_freqbaseline(cfg,freqout);


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
    
%     %% Plot TFR
%     if plot_it
%         %     cfg = [];
%         %     figure;
%         %     subplot(2,1,1);
%         %     cfg.channel      = 'OFC';
%         %     ft_singleplotTFR(cfg,freqout);
%         %     xlim([-1.5 1.5]);
%         %     subplot(2,1,2);
%         %     cfg.channel      = 'LFP';
%         %     ft_singleplotTFR(cfg,freqout);
%         %     xlim([-1.5 1.5]);
% %         
% %         cfg = [];
% %         figure;
% %         subplot(2,1,1);
% %         cfg.channel      = 'OFC';
% %         ft_singleplotTFR(cfg,freqoutbl);
% %         xlim([-0.5 2]);
% %         ylim([ 2 35]);
% %         subplot(2,1,2);
% %         cfg.channel      = 'LFP';
% %         ft_singleplotTFR(cfg,freqoutbl);
% %         xlim([-0.5 2]);
% %         ylim([ 2 35]);
%     end
%     
