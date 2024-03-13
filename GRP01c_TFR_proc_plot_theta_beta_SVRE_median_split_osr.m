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

an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr';
stat_id = 'S5t15_bhvz_nrl0_out3';
comb_method = 'mn';     % 'mn' for mean or 'md' for median
norm_wi_sbj = 1;        % 0 for TFR data, 1 for normalize within SBJ

% Plotting parameters
plot_boxes = 1;
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
if ~isempty(stat_id)
    fig_dir   = [prj_dir 'results/TFR/' an_id 'LMM/' stat_id '/median_split/'];
    eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
else
    fig_dir   = [prj_dir 'results/TFR/' an_id '/median_split/'];
end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if ~contains(an_id,'osr'), error('only run this for original sampling rate data!'); end
comb_str = ['_' comb_method];
if plot_boxes, box_str = '_boxes'; else, box_str = ''; end
if norm_wi_sbj, norm_str = '_normSBJ'; z_str = ' (z)'; else, norm_str = ''; z_str=''; end

% an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

svre_labs   = {'SV','E','pR'};
svre_names  = {'Subj. Value','Effort','Prev. Reward'};
hilo_labs   = {'High','Low'};
hilo_styles = {'-','--'};

roi_labs = {'BG','PFC'};
reds = cbrewer('seq','RdPu',5);
theta_colors = reds([2 5],:);
blues = cbrewer('seq','GnBu',5);
beta_colors = blues([2 5],:);

%% Load group model table and toss ROI and band specific outliers
if ~isempty(stat_id)
    table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all.csv'];
    fprintf('\tLoading %s...\n',table_all_fname);
    table_all = readtable(table_all_fname);

    % Toss outliers
    pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
    out_idx_all = struct;
    for f = 1:length(pow_vars)
        out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>mean(table_all.(pow_vars{f}))...
            +(st.outlier_thresh*std(table_all.(pow_vars{f})));
        if any(out_idx_all.(pow_vars{f}))
            fprintf(2,'\t%d outliers for %s in table_all\n',sum(out_idx_all.(pow_vars{f})),pow_vars{f});
        else
            fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},st.outlier_thresh);
        end
    end
end

%% Load TFR data and average within quantiles
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc_osr.mat']);
    % Check channel index
    if ~any(strcmp(tmp.tfr.label{1},{'STN','GPi'})); error('BG is not first channel!'); end
    if ~any(strcmp(tmp.tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    % Initialize data (SBJs, chan, SVRE, median_split, time)
    if s==1
        theta_avg = nan([length(SBJs) 2 3 2 numel(tmp.tfr.time)]);
        beta_avg  = nan([length(SBJs) 2 3 2 numel(tmp.tfr.time)]);
        theta_sem = nan([length(SBJs) 2 3 2 numel(tmp.tfr.time)]);
        beta_sem  = nan([length(SBJs) 2 3 2 numel(tmp.tfr.time)]);
        time_vec = tmp.tfr.time;
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
    
    % Compute theta/beta power per level
    for ch_ix = 1:2
        % Compute band-limited power
        cfg = [];
        cfg.avgoverfreq = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.frequency = squeeze(theta_lim(s,ch_ix,:))';
        theta_pow = ft_selectdata(cfg, tmp.tfr);
        cfg.frequency = squeeze(beta_lim(s,ch_ix,:))';
        beta_pow = ft_selectdata(cfg, tmp.tfr);
        if any(isnan(theta_pow.powspctrm(:))) || any(isnan(beta_pow.powspctrm(:)))
            error('nans detected!');
        end

        % Normalize power
        if norm_wi_sbj
            theta_vals = squeeze(theta_pow.powspctrm(:,ch_ix,:,:));
            theta_pow.powspctrm(:,ch_ix,:,:) = ...
                theta_pow.powspctrm(:,ch_ix,:,:)-mean(theta_vals(:))./std(theta_vals(:));
            beta_vals = squeeze(beta_pow.powspctrm(:,ch_ix,:,:));
            beta_pow.powspctrm(:,ch_ix,:,:)  = ...
                beta_pow.powspctrm(:,ch_ix,:,:)-mean(beta_vals(:))./std(beta_vals(:));
        end

        for split_ix = 1:2
            for svre_ix = 1:3
                th_idx = median_idx(svre_ix,:)==split_ix;
                b_idx = th_idx;
                if ~isempty(stat_id)
                    th_idx(out_idx_all.([roi_labs{ch_ix} '_theta'])(table_all.sbj_n==s)) = 0;
                    b_idx(out_idx_all.([roi_labs{ch_ix} '_betalo'])(table_all.sbj_n==s)) = 0;
                end
                if strcmp(comb_method,'mn')
                    % Take average and SEM of theta and beta
                    theta_avg(s,ch_ix,svre_ix,split_ix,:) = squeeze(mean(theta_pow.powspctrm(th_idx,ch_ix,1,:),1));
                    theta_sem(s,ch_ix,svre_ix,split_ix,:) = ...
                        squeeze(std(theta_pow.powspctrm(th_idx,ch_ix,1,:),[],1))./sqrt(sum(th_idx));
                    beta_avg(s,ch_ix,svre_ix,split_ix,:) = squeeze(mean(beta_pow.powspctrm(b_idx,ch_ix,1,:),1));
                    beta_sem(s,ch_ix,svre_ix,split_ix,:) = ...
                        squeeze(std(beta_pow.powspctrm(b_idx,ch_ix,1,:),[],1))./sqrt(sum(b_idx));
                elseif strcmp(comb_method,'md')
                    % Take median and MAD of theta and beta
                    theta_avg(s,ch_ix,svre_ix,split_ix,:) = squeeze(median(theta_pow.powspctrm(th_idx,ch_ix,1,:),1));
                    theta_sem(s,ch_ix,svre_ix,split_ix,:) = squeeze(mad(theta_pow.powspctrm(th_idx,ch_ix,1,:),1,1));
                    beta_avg(s,ch_ix,svre_ix,split_ix,:) = squeeze(median(beta_pow.powspctrm(b_idx,ch_ix,1,:),1));
                    beta_sem(s,ch_ix,svre_ix,split_ix,:) = squeeze(mad(beta_pow.powspctrm(b_idx,ch_ix,1,:),1,1));
                else
                    error('unknown combination method');
                end
            end
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

%% Plot SBJ TFRs
for s = 1:length(SBJs)
    %% Plot Theta median splits for each SBJ
    fig_name = [SBJs{s} '_theta_SVRE_median_split' comb_str norm_str];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        if ch_ix==1
            ch_lab = sbj_bg_roi{s};
        else
            ch_lab = sbj_pfc_roi{s};
        end
        for svre_ix = 1:3
            subplot(2,3,ch_ix*3-(3-svre_ix)); hold on;

            var_lines = gobjects([1 2]);
            for split_ix = 1:2
                % tmp = shadedErrorBar(time_vec, squeeze(theta_avg(s,ch_ix,svre_ix,split_ix,:)),...
                %     squeeze(theta_sem(s,ch_ix,svre_ix,split_ix,:)),...
                %     'lineprops',{'Color',theta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix}});
                % var_lines(split_ix) = tmp.mainLine;
                var_lines(split_ix) = plot(time_vec, squeeze(theta_avg(s,ch_ix,svre_ix,split_ix,:)), ...
                    'Color',theta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix},'LineWidth',2);
            end
            ylims = ylim;
            line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{1});

            % Axes and parameters
            title([ch_lab ' Theta Power by ' svre_names{svre_ix}], 'interpreter', 'none');
            set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
            set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
            set(gca,'YLim',ylims);
            xlabel('Time (s)');
            ylabel(['Normalized Power' z_str]);
            legend(var_lines,strcat(hilo_labs,[' ' svre_labs{svre_ix}]),'Location','best');
            set(gca,'FontSize',font_size);
        end
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end

    %% Plot Beta median splits for each SBJ
    fig_name = [SBJs{s} '_beta_SVRE_median_split' comb_str];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        if ch_ix==1
            ch_lab = sbj_bg_roi{s};
        else
            ch_lab = sbj_pfc_roi{s};
        end
        for svre_ix = 1:3
            subplot(2,3,ch_ix*3-(3-svre_ix)); hold on;

            var_lines = gobjects([1 2]);
            for split_ix = 1:2
                % tmp = shadedErrorBar(time_vec, squeeze(beta_avg(s,ch_ix,svre_ix,split_ix,:)),...
                %     squeeze(beta_sem(s,ch_ix,svre_ix,split_ix,:)),...
                %     'lineprops',{'Color',beta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix}});
                % var_lines(split_ix) = tmp.mainLine;
                var_lines(split_ix) = plot(time_vec, squeeze(beta_avg(s,ch_ix,svre_ix,split_ix,:)), ...
                    'Color',beta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix},'LineWidth',2);
            end
            ylims = ylim;
            line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{1});

            % Axes and parameters
            title([ch_lab ' Beta Power by ' svre_names{svre_ix}], 'interpreter', 'none');
            set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
            set(gca,'XTick', time_ticks);
            set(gca,'YLim',ylims);
            xlabel('Time (s)');
            ylabel(['Normalized Power' z_str]);
            legend(var_lines,strcat(hilo_labs,[' ' svre_labs{svre_ix}]),'Location','best');
            set(gca,'FontSize',font_size);
        end
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

%% Plot Group-level Theta with median splits
fig_name = ['GRP_theta_SVRE_median_split_' an_id comb_str norm_str box_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    for svre_ix = 1:3
        subplot(2,3,ch_ix*3-(3-svre_ix)); hold on;

        var_lines = gobjects([1 2]);
        for split_ix = 1:2
            var_lines(split_ix) = plot(time_vec, squeeze(theta_grp_avg(ch_ix,svre_ix,split_ix,:)), ...
                'Color',theta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix},'LineWidth',2);
            p = patch([time_vec flip(time_vec)],...
                [squeeze(theta_grp_avg(ch_ix,svre_ix,split_ix,:)+theta_grp_sem(ch_ix,svre_ix,split_ix,:))'...
                flip(squeeze(theta_grp_avg(ch_ix,svre_ix,split_ix,:)-theta_grp_sem(ch_ix,svre_ix,split_ix,:)))'],...
                theta_colors(split_ix,:),'FaceAlpha',0.1);
            p.EdgeColor = theta_colors(split_ix,:);
        end
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        if plot_boxes
            patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],...
                'k','FaceAlpha',0,'LineWidth',2,'LineStyle','--');
        end

        % Axes and parameters
        title([roi_labs{ch_ix} ' Theta Power by ' svre_names{svre_ix}], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', time_ticks);
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel(['Normalized Power' z_str]);
        legend(var_lines,strcat(hilo_labs,[' ' svre_labs{svre_ix}]),'Location','best');
        set(gca,'FontSize',font_size);
    end
end
    
% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Group-level Beta with median splits
fig_name = ['GRP_beta_SVRE_median_split_' an_id comb_str norm_str box_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for ch_ix = 1:2
    for svre_ix = 1:3
        subplot(2,3,ch_ix*3-(3-svre_ix)); hold on;

        var_lines = gobjects([1 2]);
        for split_ix = 1:2
            var_lines(split_ix) = plot(time_vec, squeeze(beta_grp_avg(ch_ix,svre_ix,split_ix,:)), ...
                'Color',beta_colors(split_ix,:),'LineStyle',hilo_styles{split_ix},'LineWidth',2);
            p = patch([time_vec flip(time_vec)],...
                [squeeze(beta_grp_avg(ch_ix,svre_ix,split_ix,:)+beta_grp_sem(ch_ix,svre_ix,split_ix,:))'...
                flip(squeeze(beta_grp_avg(ch_ix,svre_ix,split_ix,:)-beta_grp_sem(ch_ix,svre_ix,split_ix,:)))'],...
                beta_colors(split_ix,:),'FaceAlpha',0.1);
            p.EdgeColor = beta_colors(split_ix,:);
        end
        ylims = ylim;
        line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        if plot_boxes
            patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],...
                'k','FaceAlpha',0,'LineWidth',2,'LineStyle','--');
        end

        % Axes and parameters
        title([roi_labs{ch_ix} ' Beta Power by ' svre_names{svre_ix}], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', time_ticks);
        set(gca,'YLim',ylims);
        xlabel('Time (s)');
        ylabel(['Normalized Power' z_str]);
        legend(var_lines,strcat(hilo_labs,[' ' svre_labs{svre_ix}]),'Location','best');
        set(gca,'FontSize',font_size);
    end
end
    
% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
