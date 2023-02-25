%% Plots for Simon grants
clear all
close all
% clc

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_madA8t1_f2t40';%
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
outlier_thresh = 4;
use_simon_tfr = 0;
toss_same_trials = 1;

if contains(an_id,'_S') || contains(an_id,'simon')
    if contains(an_id,'A8t1')
        an_lim = [-0.8 0];
    else
        an_lim = [0.5 1.5];
    end
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end

% Analysis parameters:
theta_canon  = [4 7];
beta_canon = [12 20];
betahi_canon = [21 35];
theta_lim  = [2 8];
beta_lim   = [12 30];    
sbj_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
theta_bw = 3;
beta_bw = 4;

theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
beta_cf = [-1 -1 -1 -1; 10 17 13 12]; % PFC03, PFC04, PFC05, PFC01


% Plotting parameters
plot_psd  = 1;
font_size = 24;
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
freq_ticks = 5:5:35;
symmetric_clim = 1;
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
out_thresh_str = ['_out' num2str(outlier_thresh)];
win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');

lmm_name = [an_id win_str norm_bhv_str norm_nrl_str out_thresh_str];
lmm_fname = [prj_dir 'data/GRP/GRP_' lmm_name '_LMM_ts.mat'];
fig_dir   = [prj_dir 'results/TFR/LMM/' lmm_name '/grant_plot/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

theta_lim  = nan(length(SBJs),2,2);
beta_lim = nan(length(SBJs),2,2);
for s = 1:length(SBJs)
    for ch_ix = 1:2
        if theta_cf(ch_ix,s)==-1
            theta_lim(s,ch_ix,:) = theta_canon;
        else
            theta_lim(s,ch_ix,:) = [theta_cf(ch_ix,s)-theta_bw/2 theta_cf(ch_ix,s)+theta_bw/2];
        end
        if beta_cf(ch_ix,s)==-1
            beta_lim(s,ch_ix,:) = beta_canon;
        else
            beta_lim(s,ch_ix,:)  = [beta_cf(ch_ix,s)-beta_bw/2 beta_cf(ch_ix,s)+beta_bw/2];
        end
%         if betahi_cf(ch_ix,s)==-1
%             betahi_lim(s,ch_ix,:) = betahi_canon;
%         else
%             betahi_lim(s,ch_ix,:)  = [betahi_cf(ch_ix,s)-betahi_bw/2 betahi_cf(ch_ix,s)+betahi_bw/2];
%         end
        
        % Print final parameters
        if ch_ix==1; ch_lab = sbj_bg_roi{s}; else; ch_lab = sbj_pfc_roi{s}; end
%         fprintf('%s %s theta CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,theta_cf(ch_ix,s),theta_lim(s,ch_ix,1),theta_lim(s,ch_ix,2));
%         fprintf('%s %s beta lo CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betalo_cf(ch_ix,s),betalo_lim(s,ch_ix,1),betalo_lim(s,ch_ix,2));
%         fprintf('%s %s beta hi CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betahi_cf(ch_ix,s),betahi_lim(s,ch_ix,1),betahi_lim(s,ch_ix,2));
    end
end

%% Load Simon file details
if contains(an_id,'simon')
    DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
    DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';
    
    [numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
    % SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
    if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end
    
    FileDetails = strings(2:end,:);
end

%% Load TFR data
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    if contains(an_id,'simon')
        error('why use simon data?');
    else
        proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
        fprintf('Loading %s\n',proc_fname);
        tmp = load(proc_fname,'tfr');
    end
    
    
    % Initialize data
    if s==1
        theta_avg  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_avg   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        theta_sem  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_sem   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        time_vec = tmp.tfr.time;
        freq_vec = tmp.tfr.freq;
        tfr      = nan([length(SBJs) 2 numel(tmp.tfr.freq) numel(tmp.tfr.time)]);
        psd      = nan([length(SBJs) 2 numel(tmp.tfr.freq)]);
    end
    
    % Average across trials
    tfr(s,:,:,:) = squeeze(nanmean(tmp.tfr.powspctrm,1));
    
    % Compute PSD
    cfgs = [];
    cfgs.latency = an_lim;
    psd_tfr = ft_selectdata(cfgs,tmp.tfr);
    trl_pow = nanmean(psd_tfr.powspctrm,4);
    psd(s,:,:) = squeeze(nanmean(nanmean(psd_tfr.powspctrm,1),4));
    
    % Find PSD peaks
    for ch_ix = 1:2
        % Compute power bands
        cfg = [];
        cfg.channel = tmp.tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.frequency = squeeze(theta_lim(s,ch_ix,:))';
        theta_pow = ft_selectdata(cfg, tmp.tfr);
        theta_sem(s,ch_ix,:) = squeeze(nanstd(theta_pow.powspctrm,[],1))./sqrt(size(theta_pow.powspctrm,1));
        theta_avg(s,ch_ix,:) = squeeze(nanmean(theta_pow.powspctrm,1));
        
        cfg.frequency = beta_lim;
        beta_pow = ft_selectdata(cfg, tmp.tfr);
        beta_sem(s,ch_ix,:) = squeeze(nanstd(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
        beta_avg(s,ch_ix,:) = squeeze(nanmean(beta_pow.powspctrm,1));
    end
end

%% Load LMEs
lme_fname = [prj_dir 'data/GRP/GRP_' lmm_name '_LMM_ts.mat'];
fprintf('Loading %s\n',lme_fname);
load(lme_fname);
lme_time_vec = plt_time_vec;
lme_time_idx = nan(size(lme_time_vec));

%% Extract plotting LME data
theta_coef = nan([2 numel(lme_time_vec)]);
theta_ci   = nan([2 2 numel(lme_time_vec)]);
theta_pval = nan([2 numel(lme_time_vec)]);
theta_sig  = nan([2 numel(lme_time_vec)]);
beta_coef  = nan([2 numel(lme_time_vec)]);
beta_ci    = nan([2 2 numel(lme_time_vec)]);
beta_pval  = nan([2 numel(lme_time_vec)]);
beta_sig   = nan([2 numel(lme_time_vec)]);
for ch_ix = 1:2
    theta_coef_ix = find(strcmp(theta_lme{ch_ix,1}.CoefficientNames,'reward_prv'));
    beta_coef_ix  = find(strcmp(beta_lme{ch_ix,1}.CoefficientNames,'effortS_cur'));
    for t_ix = 1:length(plt_time_vec)
        [~, lme_time_idx(t_ix)] = min(abs(time_vec-lme_time_vec(t_ix)));
        
        theta_coef(ch_ix,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Estimate(theta_coef_ix);
        theta_pval(ch_ix,t_ix) = theta_stat{ch_ix,t_ix}.pValue(2);
        % Remove coef to get deviation for shaddedErrorBar:
        theta_ci(ch_ix,1,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Upper(theta_coef_ix)-theta_coef(ch_ix,t_ix);
        theta_ci(ch_ix,2,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Lower(theta_coef_ix)+theta_coef(ch_ix,t_ix);
        
        beta_coef(ch_ix,t_ix) = beta_lme{ch_ix,t_ix}.Coefficients.Estimate(beta_coef_ix);
        beta_pval(ch_ix,t_ix) = beta_stat{ch_ix,t_ix}.pValue(2);
        % Remove coef to get deviation for shaddedErrorBar:
        beta_ci(ch_ix,1,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Upper(theta_coef_ix)-beta_coef(ch_ix,t_ix);
        beta_ci(ch_ix,2,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Lower(theta_coef_ix)+beta_coef(ch_ix,t_ix);
    end
    
    
    % Correct for time points
    %     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [theta_sig(ch_ix,:), ~, ~, ~] = fdr_bh(theta_pval(ch_ix,:));
    [beta_sig(ch_ix,:), ~, ~, ~]  = fdr_bh(beta_pval(ch_ix,:));
end

%% Get frequency and time ticks
freq_tick_ix = nan(size(freq_ticks));
for f = 1:numel(freq_ticks)
    [~,freq_tick_ix(f)] = min(abs(freq_vec-freq_ticks(f)));
end
time_ticks = [0:0.5:2];%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
time_tick_ix = nan(size(time_ticks));
time_tick_lab = cell(size(time_ticks));
for t = 1:numel(time_ticks)
    [~,time_tick_ix(t)] = min(abs(time_vec-time_ticks(t)));
    time_tick_lab{t} = ['' num2str(time_ticks(t))];
end
if ~any(time_ticks==0); error('cant plot event with no tick at 0'); end

%% Plot SBJ TFRs
for s = 1:length(SBJs)
    %% Plot TFRs for each SBJ
    fig_name = [SBJs{s} '_TFR_theta_beta_' lmm_name];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        subplot(2,3,ch_ix*3-2); hold on;
        % Get color lims per condition
        vals = tfr(s,ch_ix,:,:);
        if symmetric_clim
            clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
            clim = [-clim clim];
        else
            clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
        end
        if ch_ix==1
            ch_lab = sbj_bg_roi{s};
            ch_color = [0 0 1];
        else
            ch_lab = sbj_pfc_roi{s};
            ch_color = [255,165,0]./255; % orange
        end
        
        % Plot TFR
        imagesc(squeeze(tfr(s,ch_ix,:,:)));
        caxis(clim);
        colorbar;
        
        % Plot Events
        line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
%         title([ch_lab '- ' an_id], 'interpreter', 'none');
        title([ch_lab ' Time-Frequency Representation'], 'interpreter', 'none');
        [~,ix1] = min(abs(time_vec-plt.plt_lim(1)));
        [~,ix2] = min(abs(time_vec-plt.plt_lim(2)));
        set(gca,'XLim',[ix1 ix2]);
        set(gca,'XTick', time_tick_ix);
        set(gca,'XTickLabels', time_ticks);
        set(gca,'YLim',[1 numel(freq_vec)]);
        set(gca,'YTick',freq_tick_ix);
        set(gca,'YTickLabels',freq_ticks);
        set(gca,'YDir','normal');
%         ylim([min(freq_vec) max(freq_vec)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        set(gca,'FontSize',font_size);
        
        %% Plot PSD
        subplot(2,3,ch_ix*3-1); hold on;
        plot(freq_vec,squeeze(psd(s,ch_ix,:)),'Color',ch_color);
%         lfp_line = plot(freq_vec,squeeze(psd(s,1,:)),'Color','b');
%         pfc_line = plot(freq_vec,squeeze(psd(s,2,:)),'Color',[255,165,0]./255); %orange
%         line(xlim,[0 0],'Color','k','LineStyle','--');
%         sig_pts = find(squeeze(sig_pow(s,ch_ix,:)));
%         scatter(freq_vec(sig_pts),zeros(size(sig_pts)),25,'k','*');
%         sig_pts = find(squeeze(sig_powq(s,ch_ix,:)));
%         scatter(freq_vec(sig_pts),zeros(size(sig_pts)),50,'r','*');
        ylims = ylim;
%         if ~isnan(theta_pk(s,ch_ix))
%             patch([thetapk_lim(s,ch_ix,1) thetapk_lim(s,ch_ix,1) ...
%                 thetapk_lim(s,ch_ix,2) thetapk_lim(s,ch_ix,2)],...
%                 [ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0.1);
%         end
%         if ~isnan(beta_pk(s,ch_ix))
%             patch([betapk_lim(s,ch_ix,1) betapk_lim(s,ch_ix,1) ...
%                 betapk_lim(s,ch_ix,2) betapk_lim(s,ch_ix,2)],...
%                 [ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0.1);
%         end
        xlabel('Frequency (Hz)');
        ylabel('Power (norm)');
%         legend([lfp_line pfc_line],{'LFP',sbj_pfc_roi{s}},'Location','best');
        title(['PSD (' num2str(an_lim(1)) '-' num2str(an_lim(2)) 's)']);
        set(gca,'FontSize',font_size);
        
        %% Plot theta and beta
        subplot(2,3,ch_ix*3); hold on;
        
        t_line = shadedErrorBar(time_vec, squeeze(theta_avg(s,ch_ix,:)),...
            squeeze(theta_sem(s,ch_ix,:)),'lineprops',{'Color','r'});
        b_line = shadedErrorBar(time_vec, squeeze(beta_avg(s,ch_ix,:)),...
            squeeze(beta_sem(s,ch_ix,:)),'lineprops',{'Color','b'});
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_lab ' Mean Power'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend([t_line.mainLine b_line.mainLine],{...
            ['theta (' num2str(theta_lim(s,ch_ix,1)) '-' num2str(theta_lim(s,ch_ix,2)) ' Hz)'], ...
            ['beta (' num2str(beta_lim(s,ch_ix,1)) '-' num2str(beta_lim(s,ch_ix,2)) ' Hz)']},'Location','best');
        set(gca,'FontSize',font_size);
        
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Average at group level
tfr_grp       = nan([2 numel(freq_vec) numel(time_vec)]);
theta_grp_avg = nan([2 numel(time_vec)]);
theta_grp_sem = nan([2 numel(time_vec)]);
beta_grp_avg  = nan([2 numel(time_vec)]);
beta_grp_sem  = nan([2 numel(time_vec)]);
for ch_ix = 1:2
    tfr_grp(ch_ix,:,:) = squeeze(nanmean(tfr(:,ch_ix,:,:),1));
    theta_grp_avg(ch_ix,:) = squeeze(mean(theta_avg(:,ch_ix,:),1));
    theta_grp_sem(ch_ix,:) = squeeze(mean(theta_sem(:,ch_ix,:),1));
    beta_grp_avg(ch_ix,:)  = squeeze(mean(beta_avg(:,ch_ix,:),1));
    beta_grp_sem(ch_ix,:)  = squeeze(mean(beta_sem(:,ch_ix,:),1));
end

%% Plot TFRs for the GROUP with coefficient time series
if plot_psd
    fig_name = ['GRP_TFR_theta_beta_' lmm_name '_withPSD'];
else
    fig_name = ['GRP_TFR_theta_beta_' lmm_name '_noPSD'];
end
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    if ch_ix==1
        ch_lab = 'Basal Ganglia';
        ch_color = [0 0 1];
    else
        ch_lab = 'Prefrontal Cortex';
        ch_color = [255,165,0]./255; % orange
    end
    
    % Plot TFR
    if plot_psd
        subplot(2,4,ch_ix*4-3); hold on;
    else
        subplot(2,3,ch_ix*3-2); hold on;
    end
    % Get color lims per condition
    vals = tfr_grp(ch_ix,:,:);
    if symmetric_clim
        clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
        clim = [-clim clim];
    else
        clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
    end
    
    imagesc(squeeze(tfr_grp(ch_ix,:,:)));
    set(gca,'YDir','normal');
    caxis(clim);
    colorbar;
    
    % Plot Events
    line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' TFR'], 'interpreter', 'none');
    % title(['LFP- ' an_id], 'interpreter', 'none');
    [~,ix1] = min(abs(time_vec-plt.plt_lim(1)));
    [~,ix2] = min(abs(time_vec-plt.plt_lim(2)));
    set(gca,'XLim',[ix1 ix2]);
    set(gca,'XTick', time_tick_ix);
    set(gca,'XTickLabels', time_tick_lab);
    set(gca,'YLim',[1 numel(freq_vec)]);
    set(gca,'YTick',freq_tick_ix);
    set(gca,'YTickLabels',freq_ticks);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    set(gca,'FontSize',font_size);
    
    % Plot PSD
    if plot_psd
        subplot(2,4,ch_ix*4-2); hold on;
        plot(freq_vec,squeeze(nanmean(psd(:,ch_ix,:),1)),'Color',ch_color);
        %     line(xlim,[0 0],'Color','k','LineStyle','--');
        xlabel('Frequency (Hz)');
        ylabel('Power (z)');
        title([ch_lab ' Reactive Power (' num2str(an_lim(1)) '-' num2str(an_lim(2)) 's)']);
        set(gca,'FontSize',font_size);
    end
    
    % Plot theta and beta
    if plot_psd
        subplot(2,4,ch_ix*4-1); hold on;
    else
        subplot(2,3,ch_ix*3-1); hold on;
    end
    
    t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(ch_ix,:)), ...
        squeeze(theta_grp_sem(ch_ix,:)),'lineprops',{'Color','r','LineWidth',2});
    b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(ch_ix,:)),...
        squeeze(beta_grp_sem(ch_ix,:)), 'lineprops',{'Color','b','LineWidth',2});
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' Power'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Power (z)');
    legend([t_line.mainLine b_line.mainLine],{...
        ['Theta (' num2str(mean(theta_lim(:,ch_ix,1))) '-' num2str(mean(theta_lim(:,ch_ix,2))) ' Hz)'],...
        ['Beta (' num2str(mean(beta_lim(:,ch_ix,1))) '-' num2str(mean(beta_lim(:,ch_ix,2))) ' Hz)']},'Location','best');
    set(gca,'FontSize',font_size);
    
    % Plot LMEs
    if plot_psd
        subplot(2,4,ch_ix*4); hold on;
    else
        subplot(2,3,ch_ix*3); hold on;
    end
%     t_line = plot(lme_time_vec, squeeze(theta_coef(ch_ix,:)),'Color','r','LineWidth',2);
%     b_line = plot(lme_time_vec, squeeze(beta_coef(ch_ix,:)),'Color','b','LineWidth',2);
    t_line = shadedErrorBar(lme_time_vec, squeeze(theta_coef(ch_ix,:)), ...
        squeeze(theta_ci(ch_ix,:,:)),'lineprops',{'Color','r','LineWidth',2});
    b_line = shadedErrorBar(lme_time_vec, squeeze(beta_coef(ch_ix,:)),...
        squeeze(beta_ci(ch_ix,:,:)), 'lineprops',{'Color','b','LineWidth',2});
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' Main Effects'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Fixed Effect Coefficient');
    legend([t_line.mainLine b_line.mainLine],{...
        ['Theta ~ previous reward'],...
        ['Beta ~ effort']},'Location','best');
    set(gca,'FontSize',font_size);
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot PFC and BG time-resolved coefficients together
fig_name = ['GRP_TFR_theta_beta_' lmm_name '_combined_coef_ts'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);

patch_alpha = 0.1;
pow_y_min = -1;
leg_line_length = [20 18]; % defaults are [30,18], so go smaller; left controls length, right doesn't seem to do anything
plot_boxes = 1;
an_lim     = [0.5 1.5];

pfc_theta_color = [178,24,43]./256;
pfc_beta_color  = [33,102,172]./256;
bg_theta_color  = [244,165,130]./256;
bg_beta_color   = [146,197,222]./256;

% Plot PFC and BG power
subplot(2,3,1); hold on;
pfc_t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(2,:)), ...
    squeeze(theta_grp_sem(2,:)),'patchSaturation',patch_alpha,'lineprops',{'Color',pfc_theta_color,'LineWidth',2});
pfc_b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(2,:)),...
    squeeze(beta_grp_sem(2,:)),'patchSaturation',patch_alpha, 'lineprops',{'Color',pfc_beta_color,'LineWidth',2});
bg_t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(1,:)), ...
    squeeze(theta_grp_sem(2,:)),'patchSaturation',patch_alpha,'lineprops',{'Color',bg_theta_color,'LineWidth',2});
bg_b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(1,:)),...
    squeeze(beta_grp_sem(2,:)),'patchSaturation',patch_alpha, 'lineprops',{'Color',bg_beta_color,'LineWidth',2});
ylims = ylim;
line([0 0],[pow_y_min max(ylims)],'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});
if plot_boxes
    patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[pow_y_min ylims(2) ylims(2) pow_y_min],'k','FaceAlpha',0,...
        'LineWidth',2,'LineStyle','--');
    fig_name = [fig_name '_boxes'];
end

% Axes and parameters
title('PFC and BG Power', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',[pow_y_min max(ylims)]);
xlabel('Time (s)');
ylabel('Power (z)');
set(gca,'FontSize',font_size);

% Plot PFC and BG power for legend
subplot(2,3,5); hold on;
pfc_t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(2,:)), ...
    squeeze(theta_grp_sem(2,:)),'patchSaturation',patch_alpha,'lineprops',{'Color',pfc_theta_color,'LineWidth',2});
pfc_b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(2,:)),...
    squeeze(beta_grp_sem(2,:)),'patchSaturation',patch_alpha, 'lineprops',{'Color',pfc_beta_color,'LineWidth',2});
bg_t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(1,:)), ...
    squeeze(theta_grp_sem(2,:)),'patchSaturation',patch_alpha,'lineprops',{'Color',bg_theta_color,'LineWidth',2});
bg_b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(1,:)),...
    squeeze(beta_grp_sem(2,:)),'patchSaturation',patch_alpha, 'lineprops',{'Color',bg_beta_color,'LineWidth',2});
ylims = ylim;
line([0 0],[pow_y_min max(ylims)],'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title('PFC and BG Power', 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',[pow_y_min max(ylims)]);
xlabel('Time (s)');
ylabel('Power (z)');
leg = legend([pfc_t_line.mainLine pfc_b_line.mainLine bg_t_line.mainLine bg_b_line.mainLine],{...
    'PFC Theta','PFC Beta','BG Theta','BG Beta'},'Location','southoutside');
%     ['PFC Theta (' num2str(mean(theta_lim(:,ch_ix,1))) '-' num2str(mean(theta_lim(:,ch_ix,2))) ' Hz)'],...
%     ['PFC Beta (' num2str(mean(beta_lim(:,ch_ix,1))) '-' num2str(mean(beta_lim(:,ch_ix,2))) ' Hz)']},'Location','best');
leg.Orientation = 'horizontal';
leg.ItemTokenSize = leg_line_length;
set(gca,'FontSize',font_size);

% Plot PFC and BG coefficient time series
subplot(2,3,2); hold on;
pfc_t_line = plot(lme_time_vec, squeeze(theta_coef(2,:)),'Color',pfc_theta_color,'LineWidth',2);
pfc_b_line = plot(lme_time_vec, squeeze(beta_coef(2,:)),'Color',pfc_beta_color,'LineWidth',2);
bg_t_line = plot(lme_time_vec, squeeze(theta_coef(1,:)),'Color',bg_theta_color,'LineWidth',2);%,'LineStyle',':');
bg_b_line = plot(lme_time_vec, squeeze(beta_coef(1,:)),'Color',bg_beta_color,'LineWidth',2);%,'LineStyle',':');
ylims = ylim;
line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});
if plot_boxes
    patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0,...
        'LineWidth',2,'LineStyle','--');
end

% Axes and parameters
title(['PFC and BG Main Effects'], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',ylims);
xlabel('Time (s)');
ylabel('Fixed Effect Coefficient');
% legend([pfc_t_line pfc_b_line bg_t_line bg_b_line],{...
%     'PFC Theta ~ previous reward','PFC Beta ~ effort', ...
%     'BG Theta ~ previous reward','BG Beta ~ effort', ...
%     },'Location','eastoutside');
set(gca,'FontSize',font_size);

% Plot again for legend
subplot(2,3,6); hold on;
pfc_t_line = plot(lme_time_vec, squeeze(theta_coef(2,:)),'Color',pfc_theta_color,'LineWidth',2);
pfc_b_line = plot(lme_time_vec, squeeze(beta_coef(2,:)),'Color',pfc_beta_color,'LineWidth',2);
bg_t_line = plot(lme_time_vec, squeeze(theta_coef(1,:)),'Color',bg_theta_color,'LineWidth',2);%,'LineStyle',':');
bg_b_line = plot(lme_time_vec, squeeze(beta_coef(1,:)),'Color',bg_beta_color,'LineWidth',2);%,'LineStyle',':');
ylims = ylim;
line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ch_lab ' Main Effects'], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',ylims);
xlabel('Time (s)');
ylabel('Fixed Effect Coefficient');
legend([pfc_t_line pfc_b_line bg_t_line bg_b_line],{...
    'PFC Theta ~ prev. REW. (p=0.04)','PFC Beta  ~ EFF. (p=0.01)', ...
    'BG Theta  ~ prev. REW. (p=0.03)','BG Beta    ~ EFF. (p=0.04)', ...
    },'Location','eastoutside');
set(gca,'FontSize',font_size);

% Save Figure
if save_fig
    fig_name = [fig_name '_legspaces'];
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot TFRs for the GROUP with SBJ PSDs
if plot_boxes; box_str = '_boxes'; else; box_str = ''; end
fig_name = ['GRP_TFR_theta_beta_' lmm_name '_withPSD_noCoef' box_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    if ch_ix==1
        ch_lab = 'Basal Ganglia';
        ch_color = [0 0 1];
    else
        ch_lab = 'Prefrontal Cortex';
        ch_color = [255,165,0]./255; % orange
    end
    
    % Plot TFR
    subplot(2,3,ch_ix*3-2); hold on;
    % Get color lims per condition
    vals = tfr_grp(ch_ix,:,:);
    if symmetric_clim
        clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
        clim = [-clim clim];
    else
        clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
    end
    
    imagesc(squeeze(tfr_grp(ch_ix,:,:)));
    set(gca,'YDir','normal');
    caxis(clim);
    colorbar;
    
    % Plot Events
    line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' TFR'], 'interpreter', 'none');
    % title(['LFP- ' an_id], 'interpreter', 'none');
    [~,ix1] = min(abs(time_vec-plt.plt_lim(1)));
    [~,ix2] = min(abs(time_vec-plt.plt_lim(2)));
    set(gca,'XLim',[ix1 ix2]);
    set(gca,'XTick', time_tick_ix);
    set(gca,'XTickLabels', time_tick_lab);
    set(gca,'YLim',[1 numel(freq_vec)]);
    set(gca,'YTick',freq_tick_ix);
    set(gca,'YTickLabels',freq_ticks);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    set(gca,'FontSize',font_size);
    
    % Plot PSDs for each SBJ
    subplot(2,3,ch_ix*3-1); hold on;
    psd_lines = gobjects(size(SBJs));
    for s = 1:length(SBJs)
        psd_lines(s) = plot(freq_vec,squeeze(nanmean(psd(s,ch_ix,:),1)),...
            'Color',sbj_colors(s,:),'LineWidth',2);
    end
    %     line(xlim,[0 0],'Color','k','LineStyle','--');
    xlabel('Frequency (Hz)');
    ylabel('Power (z)');
    if strcmp(an_id,'TFRmth_S1t2_madS8t0_f2t40') && ch_ix==2; ylim([0 12]); end
    title(['Reactive PSD (' num2str(an_lim(1)) '-' num2str(an_lim(2)) 's)']);
    legend(psd_lines,SBJs,'Location','northeast');
    set(gca,'FontSize',font_size);
    
    % Plot theta and beta
    subplot(2,3,ch_ix*3); hold on;
    
    t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(ch_ix,:)), ...
        squeeze(theta_grp_sem(ch_ix,:)),'lineprops',{'Color','r','LineWidth',2});
    b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(ch_ix,:)),...
        squeeze(beta_grp_sem(ch_ix,:)), 'lineprops',{'Color','b','LineWidth',2});
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    if plot_boxes
        patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0,...
            'LineWidth',2,'LineStyle','--');
    end
    
    % Axes and parameters
    title([ch_lab ' Power'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Power (z)');
    legend([t_line.mainLine b_line.mainLine],{...
        ['Theta (' num2str(mean(theta_lim(:,ch_ix,1))) '-' num2str(mean(theta_lim(:,ch_ix,2))) ' Hz)'],...
        ['Beta (' num2str(mean(beta_lim(:,ch_ix,1))) '-' num2str(mean(beta_lim(:,ch_ix,2))) ' Hz)']},'Location','best');
    set(gca,'FontSize',font_size);
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot outliers
% No longer relevant since I'm tossing only trials that are outliers for
% time-averaged power (from main analyses)
% pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
% plot_pows = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
% if size(outlier_ix,1)==2 && size(outlier_ix,2)==length(pow_vars) && length(size(outlier_ix))==3
%     orig_outlier_ix = outlier_ix;
%     outlier_ix = cell([length(pow_vars) length(plt_time_vec)]);
%     for t_ix = 1:length(plt_time_vec)
%         for f = 1:length(pow_vars)
%             % Track outliers per channel and frequency band
%             if contains(pow_vars{f},'PFC'), ch_ix = 2; else; ch_ix = 1; end
%             outlier_ix{f,t_ix} = find(orig_outlier_ix{ch_ix,f,t_ix});
%         end
%     end
% end
% 
% n_outliers = nan([length(pow_vars) length(plt_time_vec)]);
% for t_ix = 1:length(plt_time_vec)
%     for f = 1:length(pow_vars)
%         n_outliers(f,t_ix) = length(outlier_ix{f,t_ix});
%     end
% end
% 
% figure;
% lines = struct;
% for f = 1:length(plot_pows)
%     pow_ix = strcmp(pow_vars,plot_pows{f});
%     if contains(pow_vars{pow_ix},'BG'); ch_ix = 1; else; ch_ix = 2; end
%     if contains(pow_vars{pow_ix},'theta')
%         line_color = 'r';
%     elseif contains(pow_vars{pow_ix},'betalo')
%         line_color = 'b';
%     else
%         error('only theat and betalo now');
%     end
%     subplot(1,2,ch_ix); hold on;
%     lines.(pow_vars{pow_ix}) = plot(lme_time_vec,squeeze(n_outliers(pow_ix,:)),'Color',line_color);
% end
% subplot(1,2,1);
% title(['BG Outlier counts'], 'interpreter', 'none');
% set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
% set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
% xlabel('Time (s)');
% ylabel('# Outliers');
% legend([lines.BG_theta lines.BG_betalo],{'BG theta','BG betalo'},'Location','best');
% set(gca,'FontSize',font_size);
% subplot(1,2,2);
% title(['PFC Outlier counts'], 'interpreter', 'none');
% set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
% set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
% xlabel('Time (s)');
% ylabel('# Outliers');
% legend([lines.PFC_theta lines.PFC_betalo],{'PFC theta','PFC betalo'},'Location','best');
% set(gca,'FontSize',font_size);

