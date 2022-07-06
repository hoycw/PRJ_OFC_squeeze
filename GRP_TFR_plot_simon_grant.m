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
sbj_bg_roi   = {'GPi','STN','GPi','STN'};
man_trl_rej_ix = {[], [71 72], [], [27 28 79 80 86 87 97 98 102 103 128 139 140 148 149 150]};

% an_id = 'TFRw_S25t2_dbS25t05_fl2t40_c7';%'TFRw_S25t2_noBsln_fl1t40_c7';%'TFRw_S25t2_noBsln_fl2t40_c7';
an_id = 'simon_S';
use_simon_tfr = 1;
toss_same_trials = 1;

if contains(an_id,'_S') || contains(an_id,'simon')
    psd_win_lim = [0.5 1.5];
    peak_sign = 1;
%     an_lim = [0.5 1];
elseif contains(an_id,'_D')
    psd_win_lim = [-0.5 0];
    peak_sign = -1;
%     an_lim = [-0.5 0];
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
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

freq_ticks = 5:5:35;
symmetric_clim = 1;
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
fig_dir   = [prj_dir 'results/TFR/' an_id '/grant_plot/'];
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
        % Load data files
        proc_fname = strcat(DataStorage2,'/',FileDetails{s,1},'Stimulus_Locked.mat');
        fprintf('Loading %s\n',proc_fname);
        load(proc_fname);
        tmp.tfr = AllData.TFbl;
        bhvs{s} = AllData.exp;
        
        % Trim to plotting
        cfgs = [];
        cfgs.latency = plt.plt_lim;
        cfgs.frequency = [min(tmp.tfr.freq) 40];
        tmp.tfr = ft_selectdata(cfgs,tmp.tfr);
        
        % Remove bad trials
        if toss_same_trials
            load([sbj_dir SBJs{s} '_stim_preproc.mat'],'sbj_data');
            % Combine bad behavioral and neural trials 
            %   (empty and bad key are already tossed in Simon data)
            simon_trl_idx = 1:150;
            if strcmp(SBJs{s},'PFC04')
                simon_trl_idx(72) = [];% from whrc neural variance
            elseif strcmp(SBJs{s},'PFC05')
                simon_trl_idx(76) = [];% from whrempty
            elseif strcmp(SBJs{s},'PFC01')
                simon_trl_idx(26) = [];% from whrtrl
            end
            all_bad_ix = unique([man_trl_rej_ix{s}'; sbj_data.bhv.empty_ix; sbj_data.bhv.bad_resp_ix; sbj_data.bhv.bad_rt_ix]);
            simon_bad_ix = [];
            for t = 1:length(all_bad_ix)
                simon_bad_ix = [simon_bad_ix; find(simon_trl_idx==all_bad_ix(t))];
            end
            if isempty(simon_bad_ix)
                fprintf('%s: All bad trials already tossed!\n',SBJs{s});
            else
                fprintf(2,'%s: Removing %d trials!\n',SBJs{s},length(simon_bad_ix));
                % Remove from behavior
                bhv_fields = fieldnames(bhvs{s});
                for f = 1:length(bhv_fields)
                    if length(bhvs{s}.(bhv_fields{f}))==length(simon_trl_idx) && ~contains(bhv_fields{f},'_ix')
                        bhvs{s}.(bhv_fields{f})(simon_bad_ix) = [];
                    end
                end
                % Remove from neural
                tmp.tfr.powspctrm(simon_bad_ix,:,:,:) = [];
            end
        end
    else
%         proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
%         fprintf('Loading %s\n',proc_fname);
%         tmp = load(proc_fname,'tfr');
%         load([sbj_dir SBJs{s} '_stim_preproc.mat']);
%         bhvs{s} = sbj_data.bhv;
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
    cfgs.latency = psd_win_lim;
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
lme_fname = [prj_dir 'data/GRP/GRP_' an_id '_LME_ts.mat'];
fprintf('Loading %s\n',lme_fname);
load(lme_fname);
lme_time_vec = plt_time_vec;
lme_time_idx = nan(size(lme_time_vec));

theta_coef_ix = find(strcmp(theta_lme{1,1}.CoefficientNames,'SV_As'));
beta_coef_ix  = find(strcmp(beta_lme{1,1}.CoefficientNames,'SV_A'));
theta_coef = nan([2 numel(lme_time_vec)]);
theta_pval = nan([2 numel(lme_time_vec)]);
theta_sig  = nan([2 numel(lme_time_vec)]);
beta_coef  = nan([2 numel(lme_time_vec)]);
beta_pval  = nan([2 numel(lme_time_vec)]);
beta_sig   = nan([2 numel(lme_time_vec)]);
for ch_ix = 1:2
    for t_ix = 1:length(plt_time_vec)
        [~, lme_time_idx(t_ix)] = min(abs(time_vec-lme_time_vec(t_ix)));
        
        theta_coef(ch_ix,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.Estimate(theta_coef_ix);
        theta_pval(ch_ix,t_ix) = theta_lme{ch_ix,t_ix}.Coefficients.pValue(theta_coef_ix);
        
        beta_coef(ch_ix,t_ix) = beta_lme{ch_ix,t_ix}.Coefficients.Estimate(beta_coef_ix);
        beta_pval(ch_ix,t_ix) = beta_lme{ch_ix,t_ix}.Coefficients.pValue(beta_coef_ix);
    end
    
    % Correct for time points
    %     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [theta_sig(ch_ix,:), ~, ~, ~] = fdr_bh(theta_pval(ch_ix,:));
    [beta_sig(ch_ix,:), ~, ~, ~]  = fdr_bh(beta_pval(ch_ix,:));
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
        title(['PSD (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
        set(gca,'FontSize',16);
        
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

%% Plot TFRs for the GROUP
fig_name = ['GRP_TFR_theta_beta_' an_id];
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
    subplot(2,4,ch_ix*4-3); hold on;
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
    subplot(2,4,ch_ix*4-2); hold on;
    plot(freq_vec,squeeze(nanmean(psd(:,ch_ix,:),1)),'Color',ch_color);
%     line(xlim,[0 0],'Color','k','LineStyle','--');
    xlabel('Frequency (Hz)');
    ylabel('Power (z)');
    title([ch_lab ' Reactive Power (' num2str(psd_win_lim(1)) '-' num2str(psd_win_lim(2)) 's)']);
    set(gca,'FontSize',16);
    
    % Plot theta and beta
    subplot(2,4,ch_ix*4-1); hold on;
    
    t_line = shadedErrorBar(time_vec, squeeze(theta_grp_avg(ch_ix,:)), ...
        squeeze(theta_grp_sem(ch_ix,:)),'lineprops',{'Color','r'});
    b_line = shadedErrorBar(time_vec, squeeze(beta_grp_avg(ch_ix,:)),...
        squeeze(beta_grp_sem(ch_ix,:)), 'lineprops',{'Color','b'});
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' Power'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Power (z)');
    legend([t_line.mainLine b_line.mainLine],{...
        ['theta (' num2str(mean(theta_lim(:,ch_ix,1))) '-' num2str(mean(theta_lim(:,ch_ix,2))) ' Hz)'],...
        ['beta (' num2str(mean(beta_lim(:,ch_ix,1))) '-' num2str(mean(beta_lim(:,ch_ix,2))) ' Hz)']},'Location','best');
    set(gca,'FontSize',16);
    
    % Plot LMEs
    subplot(2,4,ch_ix*4); hold on;
    t_line = plot(lme_time_vec, squeeze(theta_coef(ch_ix,:)),'Color','r','LineWidth',2);
    b_line = plot(lme_time_vec, squeeze(beta_coef(ch_ix,:)),'Color','b','LineWidth',2);
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' Main Effects'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('Beta Coefficient');
    legend([t_line b_line],{...
        ['theta ~ previous subjective value'],...
        ['beta ~ subjective value']},'Location','best');
    set(gca,'FontSize',16);
end

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
