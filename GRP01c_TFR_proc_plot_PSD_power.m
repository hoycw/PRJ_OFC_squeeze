%% Plot TFR, baseline-corrected PSD, and theta/beta power for PFC and BG (SBJ and GRP)
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
fig_dir   = [prj_dir 'results/TFR/' an_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
if symmetric_clim; clim_str = '_sym'; else; clim_str = ''; end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load TFR data
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    
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
    
    % Compute power bands
    for ch_ix = 1:2
        cfg = [];
        cfg.channel = tmp.tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.frequency = squeeze(theta_lim(s,ch_ix,:))';
        theta_pow = ft_selectdata(cfg, tmp.tfr);
        theta_sem(s,ch_ix,:) = squeeze(nanstd(theta_pow.powspctrm,[],1))./sqrt(size(theta_pow.powspctrm,1));
        theta_avg(s,ch_ix,:) = squeeze(nanmean(theta_pow.powspctrm,1));
        
        cfg.frequency = squeeze(beta_lim(s,ch_ix,:))';
        beta_pow = ft_selectdata(cfg, tmp.tfr);
        beta_sem(s,ch_ix,:) = squeeze(nanstd(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
        beta_avg(s,ch_ix,:) = squeeze(nanmean(beta_pow.powspctrm,1));
    end
end

%% Get frequency and time ticks
freq_tick_ix = nan(size(freq_ticks));
for f = 1:numel(freq_ticks)
    [~,freq_tick_ix(f)] = min(abs(freq_vec-freq_ticks(f)));
end
time_ticks = plt.xticks;%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
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
    fig_name = [SBJs{s} '_TFR_PSD_theta_beta' clim_str];
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
        ylims = ylim;
        xlabel('Frequency (Hz)');
        ylabel('Power (norm)');
        title(['Reactive PSD (' num2str(an_lim(1)) '-' num2str(an_lim(2)) 's)']);
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

%% Plot TFRs for the GROUP with SBJ PSDs
if plot_boxes; box_str = '_boxes'; else; box_str = ''; end
fig_name = ['GRP_TFR_theta_beta_' an_id box_str clim_str];
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

