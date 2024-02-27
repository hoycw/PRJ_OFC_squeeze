%% Plot TFR and baseline-corrected PSD for PFC and BG (SBJ only) 
%   to check for 105 Hz positive control on amplifier saturation in PC+S
clear all
close all

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

an_id = 'TFRmth_S1t2_madS8t0_f80t120_osr';%'TFRmth_S1t2_madS8t0_f2t120';

% Plotting parameters
plot_boxes = 1;
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

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
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

an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

freq_ticks = min(cfg_tfr.foi):10:max(cfg_tfr.foi);
if min(cfg_tfr.foi)>105 || max(cfg_tfr.foi)<105
    error('Only run this for TFRs including 105 Hz!');
end

%% Load TFR data
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    
    % Initialize data
    if s==1
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
    psd(s,:,:) = squeeze(nanmean(nanmean(psd_tfr.powspctrm,1),4));
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
    fig_name = [SBJs{s} '_TFR_PSD_105Hz_check' clim_str];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 0.7],'Visible',fig_vis);
    
    for ch_ix = 1:2
        subplot(2,2,ch_ix*2-1); hold on;
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
        subplot(2,2,ch_ix*2); hold on;
        plot(freq_vec,squeeze(psd(s,ch_ix,:)),'Color',ch_color);
        ylims = ylim;
        xlabel('Frequency (Hz)');
        ylabel('Power (norm)');
        title(['Normalized PSD (' num2str(an_lim(1)) '-' num2str(an_lim(2)) 's)']);
        set(gca,'FontSize',font_size);
        
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

