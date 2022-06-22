%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
clc

addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_id = 'TFRw_S25t201_zbt25t05_fl2t40_varWid';%,'simon_zbt','TFRw_S25t201_zbt25t05_fl2t40'
% an_id = 'simon_D101t201_zbt25t05';%'TFRw_D101t201_zbt25t05_fl2t40_varWid';

% Analysis parameters:
theta_lim  = [4 7];
beta_lim   = [13 30];    
sbj_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
betapk_bw = 4;

% Plotting parameters
plt_id    = 'ts_S2t2_evnts_sigLine';%'ts_D1t2_evnts_sigLine';
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

%% Prep stuff
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

fig_dir   = [prj_dir 'results/TFR/' an_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

betapk_lim = nan(length(SBJs),2);
for s = 1:length(sbj_beta_pk)
    betapk_lim(s,:) = [sbj_beta_pk(s)-betapk_bw/2 sbj_beta_pk(s)+betapk_bw/2];
end

%% Load TFR data
tfr = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load data
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    
    % Initialize data
    if s==1
        theta_avg  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_avg   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        betapk_avg = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        theta_sem  = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        beta_sem   = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        betapk_sem = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        time_vec = tmp.tfr.time;
    end
    
    % Average across trials
    tfr{s} = ft_freqdescriptives([], tmp.tfr);
    
    % Compute power bands
    cfg = [];
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt  = 'no';
    cfg.frequency = theta_lim;
    theta_pow = ft_selectdata(cfg, tmp.tfr);
    theta_sem(s,:,:) = squeeze(std(theta_pow.powspctrm,[],1))./sqrt(size(theta_pow.powspctrm,1));
    theta_avg(s,:,:) = squeeze(mean(theta_pow.powspctrm,1));
    
    cfg.frequency = beta_lim;
    beta_pow = ft_selectdata(cfg, tmp.tfr);
    beta_sem(s,:,:) = squeeze(std(beta_pow.powspctrm,[],1))./sqrt(size(beta_pow.powspctrm,1));
    beta_avg(s,:,:) = squeeze(mean(beta_pow.powspctrm,1));
    
    cfg.frequency = betapk_lim(s,:);
    betapk_pow = ft_selectdata(cfg, tmp.tfr);
    betapk_sem(s,:,:) = squeeze(std(betapk_pow.powspctrm,[],1))./sqrt(size(betapk_pow.powspctrm,1));
    betapk_avg(s,:,:) = squeeze(mean(betapk_pow.powspctrm,1));
    
%     if s==1
%         time_vec = tfr{s}.time;
%         sbj_pow = tfr{s}.powspctrm;
%     else
%         sbj_pow = cat(4,sbj_pow,tfr{s}.powspctrm);
%     end
end

%% Average at group level
cfg = [];
cfg.channel = {'LFP'};
lfp_tfr_grp    = ft_freqgrandaverage(cfg,tfr{:});
lfp_theta_grp_avg  = squeeze(mean(theta_avg(:,1,:),1));
lfp_theta_grp_sem  = squeeze(std(theta_avg(:,1,:),[],1)./sqrt(size(theta_avg(:,1,:),1)));
lfp_beta_grp_avg   = squeeze(mean(beta_avg(:,1,:),1));
lfp_beta_grp_sem   = squeeze(std(beta_avg(:,1,:),[],1)./sqrt(size(beta_avg(:,1,:),1)));
lfp_betapk_grp_avg = squeeze(mean(betapk_avg(:,1,:),1));
lfp_betapk_grp_sem = squeeze(std(betapk_avg(:,1,:),[],1)./sqrt(size(betapk_avg(:,1,:),1)));

cfg.channel = {'FPC'};
fpc_tfr_grp    = ft_freqgrandaverage(cfg,tfr{strcmp(sbj_pfc_roi,'FPC')});
fpc_theta_grp_avg  = squeeze(mean(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_theta_grp_sem  = squeeze(std(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(theta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));
fpc_beta_grp_avg   = squeeze(mean(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_beta_grp_sem   = squeeze(std(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(beta_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));
fpc_betapk_grp_avg = squeeze(mean(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1));
fpc_betapk_grp_sem = squeeze(std(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),[],1)./sqrt(size(betapk_avg(strcmp(sbj_pfc_roi,'FPC'),2,:),1)));

cfg.channel = {'OFC'};
ofc_tfr_grp    = ft_freqgrandaverage(cfg,tfr{strcmp(sbj_pfc_roi,'OFC')});
ofc_theta_grp_avg  = squeeze(mean(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_theta_grp_sem  = squeeze(std(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(theta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));
ofc_beta_grp_avg   = squeeze(mean(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_beta_grp_sem   = squeeze(std(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(beta_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));
ofc_betapk_grp_avg = squeeze(mean(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1));
ofc_betapk_grp_sem = squeeze(std(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),[],1)./sqrt(size(betapk_avg(strcmp(sbj_pfc_roi,'OFC'),2,:),1)));

%% Plot SBJ TFRs
for s = 1:length(SBJs)
    %% Plot TFRs for each SBJ
    fig_name = [SBJs{s} '_TFR_theta_beta_' an_id];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:numel(tfr{s}.label)
        subplot(length(tfr{s}.label),3,ch_ix*3-2); hold on;
        % Get color lims per condition
        vals = tfr{s}.powspctrm(ch_ix,:,:);
        clim = [min(vals(:)) max(vals(:))];
        
        % Plot TFR
        imagesc(tfr{s}.time, tfr{s}.freq, squeeze(tfr{s}.powspctrm(ch_ix,:,:)));
        set(gca,'YDir','normal');
        caxis(clim);
        colorbar;
        
        % Plot Events
        line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([tfr{s}.label{ch_ix} '- ' an_id], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        ylim([min(tfr{s}.freq) max(tfr{s}.freq)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        set(gca,'FontSize',16);
        
        %% Plot theta
        subplot(length(tfr{s}.label),3,ch_ix*3-1); hold on;
        
        shadedErrorBar(time_vec, squeeze(theta_avg(s,ch_ix,:)),squeeze(theta_sem(s,ch_ix,:)));%, 'transparent',sem_alpha);
        line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([tfr{s}.label{ch_ix} ' theta (' num2str(theta_lim(1)) '-' ...
            num2str(theta_lim(2)) ' Hz): ' an_id], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        ylabel('Normalized Power');
        set(gca,'FontSize',16);
        
        %% Plot beta
        subplot(length(tfr{s}.label),3,ch_ix*3); hold on;
        
        b_line = shadedErrorBar(time_vec, squeeze(beta_avg(s,ch_ix,:)),...
            squeeze(beta_sem(s,ch_ix,:)),'lineprops',{'Color','r'});
        bpk_line = shadedErrorBar(time_vec, squeeze(betapk_avg(s,ch_ix,:)),...
            squeeze(betapk_sem(s,ch_ix,:)),'lineprops',{'Color','b'});
        line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
        title([tfr{s}.label{ch_ix} ' beta: ' an_id], 'interpreter', 'none');
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        ylabel('Normalized Power');
        legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
            ' Hz)'],['SBJ beta (' num2str(betapk_lim(s,1)) '-' num2str(betapk_lim(s,2)) ' Hz)']},'Location','best');
        set(gca,'FontSize',16);
        
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot TFRs for the GROUP
fig_name = ['GRP_TFR_theta_beta_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);

%--------------------------------------------------------------------------
%---------------------- LFP ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,3,1); hold on;
% Get color lims per condition
vals = lfp_tfr_grp.powspctrm(1,:,:);
clim = [min(vals(:)) max(vals(:))];

imagesc(lfp_tfr_grp.time, lfp_tfr_grp.freq, squeeze(lfp_tfr_grp.powspctrm(1,:,:)));
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([lfp_tfr_grp.label{1} '- ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
ylim([min(lfp_tfr_grp.freq) max(lfp_tfr_grp.freq)]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot theta
subplot(3,3,2); hold on;

shadedErrorBar(time_vec, lfp_theta_grp_avg, lfp_theta_grp_sem);%, 'transparent',sem_alpha);
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([lfp_tfr_grp.label{1} ' theta (' num2str(theta_lim(1)) '-' ...
    num2str(theta_lim(2)) ' Hz): ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
set(gca,'FontSize',16);

% Plot beta
subplot(3,3,3); hold on;

b_line = shadedErrorBar(time_vec, lfp_beta_grp_avg, lfp_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, lfp_betapk_grp_avg, lfp_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([lfp_tfr_grp.label{1} ' beta: ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(betapk_lim(s,1)) '-' num2str(betapk_lim(s,2)) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

%--------------------------------------------------------------------------
%---------------------- FPC ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,3,4); hold on;
% Get color lims per condition
vals = fpc_tfr_grp.powspctrm(1,:,:);
clim = [min(vals(:)) max(vals(:))];

imagesc(fpc_tfr_grp.time, fpc_tfr_grp.freq, squeeze(fpc_tfr_grp.powspctrm(1,:,:)));
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([fpc_tfr_grp.label{1} '- ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
ylim([min(fpc_tfr_grp.freq) max(fpc_tfr_grp.freq)]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot theta
subplot(3,3,5); hold on;

shadedErrorBar(time_vec, fpc_theta_grp_avg, fpc_theta_grp_sem);%, 'transparent',sem_alpha);
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([fpc_tfr_grp.label{1} ' theta (' num2str(theta_lim(1)) '-' ...
    num2str(theta_lim(2)) ' Hz): ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
set(gca,'FontSize',16);

% Plot beta
subplot(3,3,6); hold on;

b_line = shadedErrorBar(time_vec, fpc_beta_grp_avg, fpc_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, fpc_betapk_grp_avg, fpc_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([fpc_tfr_grp.label{1} ' beta: ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(betapk_lim(s,1)) '-' num2str(betapk_lim(s,2)) ' Hz)']},'Location','best');
set(gca,'FontSize',16);

%--------------------------------------------------------------------------
%---------------------- OFC ----------------------
%--------------------------------------------------------------------------
% Plot TFR
subplot(3,3,7); hold on;
% Get color lims per condition
vals = ofc_tfr_grp.powspctrm(1,:,:);
clim = [min(vals(:)) max(vals(:))];

imagesc(ofc_tfr_grp.time, ofc_tfr_grp.freq, squeeze(ofc_tfr_grp.powspctrm(1,:,:)));
set(gca,'YDir','normal');
caxis(clim);
colorbar;

% Plot Events
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ofc_tfr_grp.label{1} '- ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
ylim([min(ofc_tfr_grp.freq) max(ofc_tfr_grp.freq)]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'FontSize',16);

% Plot theta
subplot(3,3,8); hold on;

shadedErrorBar(time_vec, ofc_theta_grp_avg, ofc_theta_grp_sem);%, 'transparent',sem_alpha);
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ofc_tfr_grp.label{1} ' theta (' num2str(theta_lim(1)) '-' ...
    num2str(theta_lim(2)) ' Hz): ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
set(gca,'FontSize',16);

% Plot beta
subplot(3,3,9); hold on;

b_line = shadedErrorBar(time_vec, ofc_beta_grp_avg, ofc_beta_grp_sem,'lineprops',{'Color','r'});
bpk_line = shadedErrorBar(time_vec, ofc_betapk_grp_avg, ofc_betapk_grp_sem,'lineprops',{'Color','b'});
line([0 0],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
    'LineStyle',plt.evnt_styles{1});

% Axes and parameters
title([ofc_tfr_grp.label{1} ' beta: ' an_id], 'interpreter', 'none');
set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
xlabel('Time (s)');
ylabel('Normalized Power');
legend([b_line.mainLine bpk_line.mainLine],{['beta (' num2str(beta_lim(1)) '-' num2str(beta_lim(2)) ...
    ' Hz)'],['SBJ beta (' num2str(betapk_lim(s,1)) '-' num2str(betapk_lim(s,2)) ' Hz)']},'Location','best');
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
