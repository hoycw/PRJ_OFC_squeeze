%% Plot coherence and jackknife coherence for individual electrode
%   NOT intended for PFC-BG connectivity
error('plotting cohjk is just zeros because the diagonal in coherence analyses is ones; AKA it isnt inter-trial coherence');
clear all
close all
% clc

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
conn_metric = 'cohjk';%'ampcorr';
an_id = 'TFRmth_S03t2_f2t30_fourier';%'TFRmth_S1t2_f2t40_fourier';%'TFRmth_S1t2_madS8t0_f2t40';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

% Plotting parameters
toss_out  = '';
font_size = 18;
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
if contains(an_id,'f2t30')
    freq_ticks = 5:5:30;
else
    freq_ticks = 5:5:35;
end
if ~any(strcmp(conn_metric,{'coh','cohjk'}))
    error('only run for coh and cohjk!');
    symmetric_clim = 0;
end
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
% if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
% if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
% out_thresh_str = ['_out' num2str(outlier_thresh)];
% win_lim_str = [num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2))];
% lmm_name = [an_id norm_bhv_str norm_nrl_str out_thresh_str];
% lmm_fname = [prj_dir 'data/GRP/GRP_' lmm_name '_LMM_ts.mat'];
fig_dir   = [prj_dir 'results/TFR/' conn_metric '_' an_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
if contains(an_id,'_S')
    time_ticks = [0:0.5:2];%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
elseif contains(an_id,'_D')
    time_ticks = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
end


[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
if strcmp(conn_metric,'ampcorr')
    conn_metric_name = 'Amplitude Correlation (r)';
elseif strcmp(conn_metric,'PLV')
    conn_metric_name = 'Phase Locking Value';
elseif strcmp(conn_metric,'coh')
    conn_metric_name = 'Coherence';
elseif strcmp(conn_metric,'cohjk')
    conn_metric_name = 'Jackknife Coherence';
else
    error(['unknown conn_metric: ' conn_metric]);
end

%% Find outliers tossed from main time-averaged analysis
out_idx = cell([length(SBJs) 2]);
if ~isempty(toss_out)
    toss_str = ['_' toss_out];
    if ~strcmp(conn_metric,'cohjk'); error('cannot toss outliers unless jackknife'); end
    if strcmp(toss_out,'main')
        out_an_id = 'TFRmth_S1t2_madS8t0_f2t40';
        outlier_stat_id = 'S5t15_bhvz_nrlz_out4';
        outlier_thresh = 4;
        
        % Load group model table
        table_all_fname = [prj_dir 'data/GRP/GRP_' out_an_id '_' outlier_stat_id '_full_table_all.csv'];
        fprintf('\tLoading %s...\n',table_all_fname);
        table_all_avg = readtable(table_all_fname);
        
        % Identify outliers
        pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
        out_idx_all = struct;
        for f = 1:length(pow_vars)
            % Identify outliers
            out_idx_all.(pow_vars{f}) = abs(table_all_avg.(pow_vars{f}))>outlier_thresh;
            fprintf(2,'Bad trials in table_all for %s: %d\n',pow_vars{f},sum(out_idx_all.(pow_vars{f})));
        end
        
        % Combine across PFC and BG
        for s = 1:length(SBJs)
            sbj_ix = find(table_all_avg.sbj_n==s);
            out_idx{s,1} = sum([out_idx_all.PFC_theta(sbj_ix) out_idx_all.BG_theta(sbj_ix)],2);
            out_idx{s,2}  = sum([out_idx_all.PFC_betalo(sbj_ix) out_idx_all.BG_betalo(sbj_ix)],2);
        end
    end
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
        theta_avg  = nan([length(SBJs) 2 numel(tmp.conn.time)]);
        beta_avg   = nan([length(SBJs) 2 numel(tmp.conn.time)]);
%         theta_sem  = nan([length(SBJs) 2 numel(tmp.conn.time)]);
%         beta_sem   = nan([length(SBJs) 2 numel(tmp.conn.time)]);
        time_vec = tmp.conn.time;
        freq_vec = tmp.conn.freq;
        conn      = nan([length(SBJs) 2 numel(tmp.conn.freq) numel(tmp.conn.time)]);
    end
    
    for ch_ix = 1:2
        % Find frequency range
        [~,t_lo_ix] = min(abs(freq_vec-theta_lim(s,ch_ix,1)));
        [~,t_hi_ix] = min(abs(freq_vec-theta_lim(s,ch_ix,2)));
        [~,b_lo_ix] = min(abs(freq_vec-betalo_lim(s,ch_ix,1)));
        [~,b_hi_ix] = min(abs(freq_vec-betalo_lim(s,ch_ix,2)));
        
        % Check for outliers
        
        % Select data
        if strcmp(conn_metric,'coh')
            conn_data = squeeze(tmp.conn.cohspctrm(ch_ix,ch_ix,:,:));
        elseif strcmp(conn_metric,'cohjk')
            conn_data = squeeze(tmp.conn.cohspctrm(:,ch_ix,ch_ix,:,:));
            % Average excluding outliers
            if isempty(out_idx{s,ch_ix})
                if isempty(toss_out)
                    toss_str = '';
                    out_idx{s,ch_ix} = zeros([size(conn_data,1) 1]);
                else
                    out_thresh = num2str(toss_out(end));
                    error('havent figured this out yet...');
                end
            end
            conn_data = squeeze(nanmean(conn_data(~out_idx{s,ch_ix},:,:),1));
        else
            conn_data = tmp.conn.(conn_metric);
        end
        
        % Average across trials and compute connectivity within frequency bands
        conn(s,ch_ix,:,:) = conn_data;
        theta_avg(s,ch_ix,:) = squeeze(nanmean(conn_data(t_lo_ix:t_hi_ix,:),1));
        beta_avg(s,ch_ix,:)  = squeeze(nanmean(conn_data(b_lo_ix:b_hi_ix,:),1));
    end
end

%% Get frequency and time ticks
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

%% Plot SBJ connectivity TFRs
for s = 1:length(SBJs)
    for ch_ix = 1:2
        %% Plot connectivity TFRs for each SBJ
        fig_name = [SBJs{s} '_TFR_theta_beta_' an_id '_' conn_metric '_' ch_labs{s}{ch_ix} toss_str];
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 1 1],'Visible',fig_vis);
        
        subplot(1,2,1); hold on;
        % Get color lims per condition
        vals = conn(s,ch_ix,:,:);
        clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
        
        % Plot TFR
        imagesc(squeeze(conn(s,ch_ix,:,:)));
        caxis(clim);
        colorbar;
        
        % Plot Events
        line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_labs{s}{ch_ix} ': ' conn_metric_name], 'interpreter', 'none');
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
        
        t_line = plot(time_vec, squeeze(theta_avg(s,ch_ix,:)),'Color','r','LineWidth',3);
        b_line = plot(time_vec, squeeze(beta_avg(s,ch_ix,:)),'Color','b','LineWidth',3);
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_labs{s}{ch_ix} ' Mean'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel(conn_metric_name);
        legend([t_line b_line],{...
            ['theta (' num2str(theta_lim(s,ch_ix,1)) '-' num2str(theta_lim(s,ch_ix,2)) ' Hz)'], ...
            ['beta (' num2str(betalo_lim(s,ch_ix,1)) '-' num2str(betalo_lim(s,ch_ix,2)) ' Hz)']},'Location','best');
        set(gca,'FontSize',font_size);
        
        
        % Save Figure
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
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
fig_name = ['GRP_TFR_theta_beta_' conn_metric '_' an_id toss_str];
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

