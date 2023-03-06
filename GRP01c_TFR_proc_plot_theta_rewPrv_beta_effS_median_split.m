%% Plot TFR, baseline-corrected PSD, and theta/beta power for PFC and BG (SBJ and GRP)
clear all
close all
% clc

addpath('/Users/colinhoy/Code/Apps/cbrewer/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_madA8t1_f2t40';%

reds = cbrewer('seq','RdPu',5);
rew_colors = reds([2 5],:);
blues = cbrewer('seq','GnBu',5);
eff_colors = blues([2 5],:);

% Plotting parameters
plot_boxes = 0;
plot_psd   = 1;
symmetric_clim = 1;
font_size  = 24;
save_fig   = 1;
fig_ftype  = 'png';
fig_vis    = 'on';

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
[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
beta_lim   = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
% betahi_lim = fn_compute_freq_lim(SBJs,betahi_cf,'betahi');

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
freq_ticks = 5:5:35;
if contains(an_id,'A8t1')
    plt_id = 'ts_S8t2_evnts_sigLine';
elseif contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
fig_dir   = [prj_dir 'results/TFR/' an_id '/median_split/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
if symmetric_clim; clim_str = '_sym'; else; clim_str = ''; end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load TFR data and average within quantiles
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    
    % Initialize data
    if s==1
        theta_avg = nan([length(SBJs) 2 2 numel(tmp.tfr.time)]);
        beta_avg  = nan([length(SBJs) 2 2 numel(tmp.tfr.time)]);
        theta_sem = nan([length(SBJs) 2 2 numel(tmp.tfr.time)]);
        beta_sem  = nan([length(SBJs) 2 2 numel(tmp.tfr.time)]);
        time_vec = tmp.tfr.time;
    end
    
    % Get stake and effort splitd (0 = low, 1 = high)
    stake_idx = sbj_data.bhv.stake>=median(sbj_data.bhv.stake);
    eff_idx   = sbj_data.bhv.effort>=median(sbj_data.bhv.effort);
    
    % Compute theta power per stake level
    for ch_ix = 1:2
        cfg = [];
        cfg.channel = tmp.tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.frequency = squeeze(theta_lim(s,ch_ix,:))';
        for stk_ix = 1:2
            cfg.trials = find(stake_idx==stk_ix-1);
            theta_pow = ft_selectdata(cfg, tmp.tfr);
            theta_sem(s,ch_ix,stk_ix,:) = squeeze(nanstd(theta_pow.powspctrm,[],1))./sqrt(size(theta_pow.powspctrm,1));
            theta_avg(s,ch_ix,stk_ix,:) = squeeze(nanmean(theta_pow.powspctrm,1));
        end
        
        cfg.frequency = squeeze(theta_lim(s,ch_ix,:))';
        for eff_ix = 1:2
            cfg.trials = find(eff_idx==eff_ix-1);
            beta_pow = ft_selectdata(cfg, tmp.tfr);
            beta_sem(s,ch_ix,eff_ix,:) = squeeze(nanstd(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
            beta_avg(s,ch_ix,eff_ix,:) = squeeze(nanmean(beta_pow.powspctrm,1));
        end
    end
end

%% Get frequency and time ticks
time_ticks = plt.xticks;%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
time_tick_ix = nan(size(time_ticks));
time_tick_lab = cell(size(time_ticks));
for t = 1:numel(time_ticks)
    [~,time_tick_ix(t)] = min(abs(time_vec-time_ticks(t)));
    time_tick_lab{t} = ['' num2str(time_ticks(t))];
end
if ~any(time_ticks==0); error('cant plot event with no tick at 0'); end

rew_lab = {'Low Prev. Reward','High Prev. Reward'};
eff_lab = {'Low Effort','High Effort'};

%% Plot SBJ TFRs
for s = 1:length(SBJs)
    %% Plot TFRs for each SBJ
    fig_name = [SBJs{s} '_TFR_median_theta_rewPrv_beta_effS'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        %% Plot theta per previous reward level
        subplot(2,2,ch_ix*2-1); hold on;
        if ch_ix==1
            ch_lab = sbj_bg_roi{s};
        else
            ch_lab = sbj_pfc_roi{s};
        end
        
        r_lines = gobjects([1 2]);
        for stk_ix = 1:2
%             tmp = shadedErrorBar(time_vec, squeeze(theta_avg(s,ch_ix,stk_ix,:)),...
%                 squeeze(theta_sem(s,ch_ix,stk_ix,:)),'lineprops',{'Color',rew_colors(stk_ix,:)});
%             r_lines(stk_ix) = tmp.mainLine;
            r_lines(stk_ix) = plot(time_vec, squeeze(theta_avg(s,ch_ix,stk_ix,:)), ...
                'Color',rew_colors(stk_ix,:),'LineWidth',2);
        end
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_lab ' Theta Power by Prev Reward'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend(r_lines,rew_lab,'Location','best');
        set(gca,'FontSize',font_size);
        
        %% Plot beta per effort level
        subplot(2,2,ch_ix*2); hold on;
        e_lines = gobjects([1 2]);
        for eff_ix = 1:2
%             tmp = shadedErrorBar(time_vec, squeeze(beta_avg(s,ch_ix,eff_ix,:)),...
%                 squeeze(beta_sem(s,ch_ix,eff_ix,:)),'lineprops',{'Color',eff_colors(eff_ix,:)});
%             e_lines(eff_ix) = tmp.mainLine;
            e_lines(eff_ix) = plot(time_vec, squeeze(beta_avg(s,ch_ix,eff_ix,:)), ...
                'Color',eff_colors(eff_ix,:),'LineWidth',2);
        end
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([ch_lab ' Low Beta Power by Effort'], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend(e_lines,eff_lab,'Location','best');
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
theta_grp_avg = squeeze(mean(theta_avg,1));
theta_grp_sem = squeeze(mean(theta_sem,1));
beta_grp_avg  = squeeze(mean(beta_avg,1));
beta_grp_sem  = squeeze(mean(beta_sem,1));

%% Plot TFRs for the GROUP with SBJ PSDs
if plot_boxes; box_str = '_boxes'; else; box_str = ''; end
fig_name = ['GRP_TFR_median_theta_rewPrv_beta_effS_' an_id box_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    if ch_ix==1
        ch_lab = 'Basal Ganglia';
    else
        ch_lab = 'Prefrontal Cortex';
    end
    
    % Plot theta by previous reward
    subplot(2,2,ch_ix*2-1); hold on;
    
    r_lines = gobjects([1 2]);
    for stk_ix = 1:2
        r_lines(stk_ix) = plot(time_vec, squeeze(theta_grp_avg(ch_ix,stk_ix,:)), ...
            'Color',rew_colors(stk_ix,:),'LineWidth',2);
        p = patch([time_vec flip(time_vec)],[squeeze(theta_grp_avg(ch_ix,stk_ix,:)+theta_grp_sem(ch_ix,stk_ix,:))'...
            flip(squeeze(theta_grp_avg(ch_ix,stk_ix,:)-theta_grp_sem(ch_ix,stk_ix,:)))'],...
            rew_colors(stk_ix,:),'FaceAlpha',0.1);
        p.EdgeColor = rew_colors(stk_ix,:);
    end
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    if plot_boxes
        patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0,...
            'LineWidth',2,'LineStyle','--');
    end
    
    % Axes and parameters
    title([ch_lab ' Theta Power by Previous Reward'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Power (z)');
    legend(r_lines,rew_lab,'Location','best');
    set(gca,'FontSize',font_size);
    
    % Plot beta by effort
    subplot(2,2,ch_ix*2); hold on;
    
    e_lines = gobjects([1 2]);
    for eff_ix = 1:2
        e_lines(eff_ix) = plot(time_vec, squeeze(beta_grp_avg(ch_ix,eff_ix,:)), ...
            'Color',eff_colors(eff_ix,:),'LineWidth',2);
        p = patch([time_vec flip(time_vec)],[squeeze(beta_grp_avg(ch_ix,eff_ix,:)+beta_grp_sem(ch_ix,eff_ix,:))' ...
            flip(squeeze(beta_grp_avg(ch_ix,eff_ix,:)-beta_grp_sem(ch_ix,eff_ix,:)))'],...
            eff_colors(eff_ix,:),'FaceAlpha',0.1);
        p.EdgeColor = eff_colors(eff_ix,:);
    end
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    if plot_boxes
        patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0,...
            'LineWidth',2,'LineStyle','--');
    end
    
    % Axes and parameters
    title([ch_lab ' Low Beta Power by Effort'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Power (z)');
    legend(e_lines,eff_lab,'Location','best');
    set(gca,'FontSize',font_size);
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

