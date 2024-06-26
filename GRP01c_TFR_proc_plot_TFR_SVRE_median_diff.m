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

an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr';%'TFRmth_S1t2_madA8t1_f2t40';%

% Plotting parameters
plot_boxes = 1;
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
fig_dir   = [prj_dir 'results/TFR/' an_id '/median_diff/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
if ~contains(an_id,'osr'), error('only run this for original sampling rate data!'); end
if plot_boxes, box_str = '_boxes'; else, box_str = ''; end
if symmetric_clim; clim_str = '_sym'; else; clim_str = ''; end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

svre_labs   = {'SV','E','pR'};
svre_names  = {'Subj. Value','Effort','Prev. Reward'};

%% Load TFR data
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc_osr.mat']);
    
    % Initialize data
    if s==1
        time_vec = tmp.tfr.time;
        freq_vec = tmp.tfr.freq;
        tfr      = nan([length(SBJs) 2 3 2 numel(tmp.tfr.freq) numel(tmp.tfr.time)]);
    end
    
    % Get SV, reward, and effort splits
    % (1 = low, 2 = high; exclude middle level = 0 for R/E)
    median_idx = zeros([length(svre_labs) length(sbj_data.bhv.SV)]);
    median_idx(1,sbj_data.bhv.SV<=median(sbj_data.bhv.SV)) = 1;
    median_idx(1,sbj_data.bhv.SV>median(sbj_data.bhv.SV)) = 2;
    median_idx(2,sbj_data.bhv.effort<0.4) = 1;
    median_idx(2,sbj_data.bhv.effort>0.5) = 2;
    median_idx(3,sbj_data.bhv.stake_prv<6) = 1;
    median_idx(3,sbj_data.bhv.stake_prv>8) = 2;

    % Average across trials
    for split_ix = 1:2
        for svre_ix = 1:3
            tfr(s,:,svre_ix,split_ix,:,:) = squeeze(nanmean(tmp.tfr.powspctrm(median_idx(svre_ix,:)==split_ix,:,:,:),1));
        end
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

%% Average at group level
% Take difference across median splits
tfr_diff = squeeze(diff(tfr,1,4));

tfr_grp       = nan([2 3 numel(freq_vec) numel(time_vec)]);
tfr_diff_grp = squeeze(nanmean(tfr_diff,1));

%% Plot TFRs for the GROUP
if plot_boxes; box_str = '_boxes'; else; box_str = ''; end
fig_name = ['GRP_TFR_theta_beta_' an_id box_str clim_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    if ch_ix==1
        ch_lab = 'Basal Ganglia';
    else
        ch_lab = 'Prefrontal Cortex';
    end

    for svre_ix = 1:3
        % Plot TFR
        subplot(2,3,ch_ix*3-(3-svre_ix)); hold on;
        % Get color lims per condition
        vals = tfr_diff_grp(ch_ix,svre_ix,:,:);
        if symmetric_clim
            clim = max(abs([prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))]));
            clim = [-clim clim];
        else
            clim = [prctile(vals(:),plt.clim_perc(1)) prctile(vals(:),plt.clim_perc(2))];
        end

        imagesc(squeeze(tfr_diff_grp(ch_ix,svre_ix,:,:)));
        set(gca,'YDir','normal');
        caxis(clim);
        cbar = colorbar;
        cbar.Label.String = ['Hi-Lo: ' svre_names{svre_ix}];

        % Plot Events
        line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});

        % Plot boxes
        if plot_boxes
            [~,ix1] = min(abs(time_vec-an_lim(1)));
            [~,ix2] = min(abs(time_vec-an_lim(2)));
            patch([ix1 ix1 ix2 ix2],...
                [min(theta_lim(:)) max(theta_lim(:)) max(theta_lim(:)) min(theta_lim(:))],...
                'k','EdgeColor',[1 1 1],'FaceAlpha',0,'LineWidth',3,'LineStyle','--');
            patch([ix1 ix1 ix2 ix2],...
                [min(beta_lim(:)) max(beta_lim(:)) max(beta_lim(:)) min(beta_lim(:))],...
                'k','EdgeColor',[1 1 1],'FaceAlpha',0,'LineWidth',3,'LineStyle','--');
        end

        % Axes and parameters
        title([ch_lab ' TFR: ' svre_names{svre_ix}], 'interpreter', 'none');
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
    end
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

