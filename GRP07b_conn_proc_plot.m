%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
conn_metric = 'ampcorr';
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_madS8t0_f2t40';
freq_ch = 'PFC';    % which channel's SBJ-specific frequency band should be used? 'PFC' or 'BG'

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

% Plotting parameters
font_size = 18;
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
freq_ticks = 5:5:35;
if strcmp(conn_metric,'PLV')
    symmetric_clim = 0;
else
    symmetric_clim = 1;
end
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
    time_ticks = [0:0.5:2];%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
    time_ticks = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
else
    error('couldnt pick plt_id based on an_id');
end
% if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
% if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
% out_thresh_str = ['_out' num2str(outlier_thresh)];
% win_lim_str = [num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2))];
% lmm_name = [an_id norm_bhv_str norm_nrl_str out_thresh_str];
% lmm_fname = [prj_dir 'data/GRP/GRP_' lmm_name '_LMM_ts.mat'];
fig_dir   = [prj_dir 'results/TFR/' conn_metric '_' an_id '/grant_plot/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
if strcmp(freq_ch,'BG')
    pk_frq_ch_ix = 1;
elseif strcmp(freq_ch,'PFC')
    pk_frq_ch_ix = 2;
else
    error('unknown peak frequency channel');
end
if strcmp(conn_metric,'ampcorr')
    conn_metric_name = 'Amplitude Correlation (r)';
elseif strcmp(conn_metric,'PLV')
    conn_metric_name = 'Phase Locking Value';
else
    error(['unknown conn_metric: ' conn_metric]);
end

%% Load connectivity data
ch_labs = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load connectivity
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '_' conn_metric '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'conn');
    ch_labs{s} = tmp.conn.label;
    
    % Initialize data
    if s==1
        theta_avg  = nan([length(SBJs) numel(tmp.conn.time)]);
        beta_avg   = nan([length(SBJs) numel(tmp.conn.time)]);
%         theta_sem  = nan([length(SBJs) numel(tmp.conn.time)]);
%         beta_sem   = nan([length(SBJs) numel(tmp.conn.time)]);
        time_vec = tmp.conn.time;
        freq_vec = tmp.conn.freq;
        conn      = nan([length(SBJs) numel(tmp.conn.freq) numel(tmp.conn.time)]);
    end
    
    % Average across trials
    conn(s,:,:) = tmp.conn.(conn_metric);
    
    % Compute connectivity within frequency bands
    [~,lo_ix] = min(abs(freq_vec-theta_lim(s,pk_frq_ch_ix,1)));
    [~,hi_ix] = min(abs(freq_vec-theta_lim(s,pk_frq_ch_ix,2)));
%     theta_sem(s,:) = squeeze(nanstd(tmp.conn.(conn_metric)(lo_ix:hi_ix,:),[],1)) ...
%               ./sqrt(size(tmp.conn.(conn_metric)(lo_ix:hi_ix,:),1));
    theta_avg(s,:) = squeeze(nanmean(tmp.conn.(conn_metric)(lo_ix:hi_ix,:),1));
    
    [~,lo_ix] = min(abs(freq_vec-betalo_lim(s,pk_frq_ch_ix,1)));
    [~,hi_ix] = min(abs(freq_vec-betalo_lim(s,pk_frq_ch_ix,2)));
%     beta_sem(s,:) = squeeze(nanstd(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
    beta_avg(s,:) = squeeze(nanmean(tmp.conn.(conn_metric)(lo_ix:hi_ix,:),1));
end

%% Plot SBJ connectivity TFRs
% Get frequency and time ticks
freq_tick_ix = nan(size(freq_ticks));
for f = 1:numel(freq_ticks)
    [~,freq_tick_ix(f)] = min(abs(freq_vec-freq_ticks(f)));
end
time_tick_ix = nan(size(time_ticks));
time_tick_lab = cell(size(time_ticks));
for t = 1:numel(time_ticks)
    [~,time_tick_ix(t)] = min(abs(time_vec-time_ticks(t)));
    time_tick_lab{t} = ['' num2str(time_ticks(t))];
end
if ~any(time_ticks==0); error('cant plot event with no tick at 0'); end

for s = 1:length(SBJs)
    %% Plot connectivity TFRs for each SBJ
    fig_name = [SBJs{s} '_TFR_theta_beta_' an_id '_' conn_metric];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    subplot(1,2,1); hold on;
    % Get color lims per condition
    vals = conn(s,:,:);
    if symmetric_clim
        clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
        clim = [-clim clim];
    else
        clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
    end
    
    % Plot TFR
    imagesc(squeeze(conn(s,:,:)));
    caxis(clim);
    colorbar;
    
    % Plot Events
    line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_labs{s} ' Connectivity: ' conn_metric_name], 'interpreter', 'none');
    [~,ix1] = min(abs(time_vec-plt.plt_lim(1)));
    [~,ix2] = min(abs(time_vec-plt.plt_lim(2)));
    set(gca,'XLim',[ix1 ix2]);
%     set(gca,'XLim',[1 numel(time_vec)]);
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
    
    %% Plot connectivity averaged within theta and beta
    subplot(1,2,2); hold on;
    
    t_line = plot(time_vec, squeeze(theta_avg(s,:)),'Color','r','LineWidth',3);
    b_line = plot(time_vec, squeeze(beta_avg(s,:)),'Color','b','LineWidth',3);
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_labs{s} ' Mean Connectivity'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel(conn_metric_name);
    legend([t_line b_line],{...
        ['theta (' num2str(theta_lim(s,pk_frq_ch_ix,1)) '-' num2str(theta_lim(s,pk_frq_ch_ix,2)) ' Hz)'], ...
        ['beta (' num2str(betalo_lim(s,pk_frq_ch_ix,1)) '-' num2str(betalo_lim(s,pk_frq_ch_ix,2)) ' Hz)']},'Location','best');
    set(gca,'FontSize',font_size);
    
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Average at group level
conn_grp = squeeze(mean(conn,1));

% Summary stats for theta
theta_grp_avg = squeeze(mean(theta_avg,1));
theta_grp_sem = squeeze(std(theta_avg,1))./sqrt(length(SBJs));
% Summary stats for low beta
beta_grp_avg  = squeeze(mean(beta_avg,1));
beta_grp_sem  = squeeze(std(beta_avg,1))./sqrt(length(SBJs));

%% Plot TFRs for the GROUP
fig_name = ['GRP_TFR_theta_beta_' conn_metric '_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
ch_lab = 'PFC-BG';

% Get color lims per condition
vals = conn_grp;
if symmetric_clim
    clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
    clim = [-clim clim];
else
    clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
end

% Plot connectivity TFR
subplot(1,2,1); hold on;
imagesc(conn_grp);
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
    'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ch_lab ' Connectivity: ' conn_metric_name], 'interpreter', 'none');
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

% Plot theta and beta
subplot(1,2,2); hold on;
t_line = shadedErrorBar(time_vec, theta_grp_avg, ...
    theta_grp_sem,'lineprops',{'Color','r'});
b_line = shadedErrorBar(time_vec, beta_grp_avg,...
    beta_grp_sem, 'lineprops',{'Color','b'});
ylims = ylim;
line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ch_lab ' Band-Specific ' conn_metric_name], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',ylims);
xlabel('Time (s)');
ylabel(conn_metric_name);
legend([t_line.mainLine b_line.mainLine],{...
    ['Theta (' num2str(mean(theta_lim(:,pk_frq_ch_ix,1))) '-' num2str(mean(theta_lim(:,pk_frq_ch_ix,2))) ' Hz)'],...
    ['Beta (' num2str(mean(betalo_lim(:,pk_frq_ch_ix,1))) '-' num2str(mean(betalo_lim(:,pk_frq_ch_ix,2))) ' Hz)']},'Location','best');
set(gca,'FontSize',font_size);


% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

