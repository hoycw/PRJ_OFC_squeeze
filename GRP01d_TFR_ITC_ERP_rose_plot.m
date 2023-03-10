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

an_id = 'TFRmth_S1t2_f2t40_fourier';

% Plotting parameters
plot_boxes = 0;
plot_psd   = 1;
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
fig_dir   = [prj_dir 'results/TFR/' an_id '/ITPC_rose/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if ~strcmp(cfg_tfr.output,'fourier'); error('need to run on fourier output for pahse analyses!'); end
plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load TFR data
theta_ang = cell([length(SBJs) 2]);
beta_ang  = cell([length(SBJs) 2]);
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    preproc_fname = [sbj_dir SBJs{s} '_stim_preproc.mat']; % This is [-3 10] to cover all possible times
    fprintf('Loading %s\n',preproc_fname);
    load(preproc_fname,'sbj_data');

    % Initialize data
    if s==1
        time_vec  = tmp.tfr.time;
        freq_vec  = tmp.tfr.freq;
        erps      = nan([length(SBJs) 2 numel(tmp.tfr.time)]);
        itpc      = nan([length(SBJs) 2 numel(tmp.tfr.freq) numel(tmp.tfr.time)]);
    end
    
    % Average ERP
    if strcmp(an.event_type,'D') && strcmp(an.bsln_evnt,'S')
        error('no time series re-aligbnment to decision implemented yet');
    else
        cfg_erp = [];
        cfg_erp.hpfilter       = 'yes';
        cfg_erp.hpfreq         = 1;
        % cfgpp.hpfiltord      = 4; % Leaving blank causes instability error, 1 or 2 works
        cfg_erp.lpfilter       = 'yes';
        cfg_erp.lpfreq         = 20;
        cfg_erp.demean         = 'yes';
        cfg_erp.baselinewindow = [-0.2 0];
        erp_sbj = ft_preprocessing(cfg_erp,sbj_data.ts);

        cfg_avg = [];
        cfg_avg.latency = an.trial_lim_s;
        cfg_avg.avgoverrpt = 'yes';
        erp_mn = ft_selectdata(cfg_avg,erp_sbj);
        erps(s,:,:) = erp_mn.trial{1};
    end
    
    % Angle Extraction for polar plot
    cfg_ang = [];
    cfg_ang.avgoverfreq = 'yes';
    cfg_ang.avgovertime = 'yes';
    cfg_ang.latency     = an_lim;
    for ch_ix = 1:2
        % Compute mean phase angle in T-F window
        cfg_ang.channel = tmp.tfr.label(ch_ix);
        cfg_ang.frequency = squeeze(theta_lim(s,ch_ix,:))';
        complex = ft_selectdata(cfg_ang,tmp.tfr);
        theta_ang{s,ch_ix} = squeeze(angle(complex.fourierspctrm));
        
        cfg_ang.frequency = squeeze(beta_lim(s,ch_ix,:))';
        complex = ft_selectdata(cfg_ang,tmp.tfr);
        beta_ang{s,ch_ix} = squeeze(angle(complex.fourierspctrm));
    end
    
    % Compute ITPC
    F = tmp.tfr.fourierspctrm;
    itpc_sbj = F./abs(F);       % Normalize to unit circle
    itpc_sbj = sum(itpc_sbj,1);     % Sum phase angles
    itpc(s,:,:,:) = abs(itpc_sbj)/size(tmp.tfr.fourierspctrm,1);     % Get mean of angles for consistency
end

%% Compute band-limited ITPC
freq_idx = freq_vec>=theta_lim(s,ch_ix,1) & freq_vec<=theta_lim(s,ch_ix,2);
theta_itpc = squeeze(mean(itpc(:,:,freq_idx,:),3));

freq_idx = freq_vec>=beta_lim(s,ch_ix,1) & freq_vec<=beta_lim(s,ch_ix,2);
beta_itpc = squeeze(mean(itpc(:,:,freq_idx,:),3));

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
    fig_name = [SBJs{s} '_TFR_ITC_ERP_rose_theta_beta'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    
    for ch_ix = 1:2
        subplot(2,4,[ch_ix*4-3 ch_ix*4-2]); hold on;
        % Get color lims per condition
        if ch_ix==1
            ch_lab = sbj_bg_roi{s};
            ch_color = [0 0 1];
        else
            ch_lab = sbj_pfc_roi{s};
            ch_color = [255,165,0]./255; % orange
        end
        
        % Compute mean ITPC within window
        time_idx = time_vec>=an_lim(1) & time_vec<=an_lim(2);
        theta_itpc_mean = mean(theta_itpc(s,ch_ix,time_idx),3);
        beta_itpc_mean  = mean(beta_itpc(s,ch_ix,time_idx),3);
        
        % Plot TFR
        imagesc(squeeze(itpc(s,ch_ix,:,:)));
        caxis([0 max(itpc(:))]);
        colorbar;
        
        % Plot Events
        line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{1});
        
        % Axes and parameters
%         title([ch_lab '- ' an_id], 'interpreter', 'none');
        title([ch_lab ' ITPC'], 'interpreter', 'none');
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
        
        % Plot ERP Mean
        yyaxis right;
        plot(squeeze(erps(s,ch_ix,:)),'Color','k','LineWidth',2,'LineStyle','-');
        ylabel('ERP Amplitude');
        
        %% Plot Theta Rose plot of phase angles
        subplot(2,4,ch_ix*4-1);
        polarhistogram(theta_ang{s,ch_ix},[-pi:pi/5:pi],'Normalization','probability');
        title(['theta ITPC avg(' num2str(an_lim(1),'%.1f') '-' ...
            num2str(an_lim(2),'%.1f') 's)=' num2str(theta_itpc_mean,'%.2f')]);
        set(gca,'FontSize',font_size);
        
        %% Plot theta and beta
        subplot(2,4,ch_ix*4);
        polarhistogram(beta_ang{s,ch_ix},[-pi:pi/5:pi],'Normalization','probability');
        title(['beta ITPC avg(' num2str(an_lim(1),'%.1f') '-' ...
            num2str(an_lim(2),'%.1f') 's)=' num2str(beta_itpc_mean,'%.2f')]);
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
itpc_grp       = squeeze(nanmean(itpc,1));
theta_itpc_grp = squeeze(mean(theta_itpc,1));
theta_itpc_grp_sem = squeeze(std(theta_itpc,[],1))./sqrt(size(theta_itpc,1));
beta_itpc_grp  = squeeze(mean(beta_itpc,1));
beta_itpc_grp_sem = squeeze(std(beta_itpc,[],1))./sqrt(size(beta_itpc,1));
erp_grp = squeeze(mean(erps,1));
erp_grp_sem = squeeze(std(erps,[],1))./sqrt(size(erps,1));

%% Plot TFRs for the GROUP with SBJ PSDs
if plot_boxes; box_str = '_boxes'; else; box_str = ''; end
fig_name = ['GRP_TFR_ITPC_theta_beta_' an_id box_str];
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
    subplot(2,2,ch_ix*2-1); hold on;
    
    imagesc(squeeze(itpc_grp(ch_ix,:,:)));
    set(gca,'YDir','normal');
    caxis([0 max(itpc_grp(:))]);
    colorbar;
    
    % Plot Events
    line([time_tick_ix(time_ticks==0) time_tick_ix(time_ticks==0)],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,'LineStyle',plt.evnt_styles{1});
    
    % Axes and parameters
    title([ch_lab ' ITPC'], 'interpreter', 'none');
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
    
    % Plot ERP
    yyaxis right;
    shadedErrorBar(1:size(itpc_grp,3),squeeze(erp_grp(ch_ix,:)), squeeze(erp_grp_sem(ch_ix,:)),...
        'lineprops',{'Color','k','LineWidth',2});
    ylabel('ERP Amplitude');
    
    %% Plot theta and beta ITPC
    subplot(2,2,ch_ix*2); hold on;
    
    t_line = shadedErrorBar(time_vec, squeeze(theta_itpc_grp(ch_ix,:)), ...
        squeeze(theta_itpc_grp_sem(ch_ix,:)),'lineprops',{'Color','r','LineWidth',2});
    b_line = shadedErrorBar(time_vec, squeeze(beta_itpc_grp(ch_ix,:)),...
        squeeze(beta_itpc_grp_sem(ch_ix,:)), 'lineprops',{'Color','b','LineWidth',2});
    ylims = ylim;
    line([0 0],ylims,'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    if plot_boxes
        patch([an_lim(1) an_lim(1) an_lim(2) an_lim(2)],[ylims(1) ylims(2) ylims(2) ylims(1)],'k','FaceAlpha',0,...
            'LineWidth',2,'LineStyle','--');
    end
    
    % Axes and parameters
    title([ch_lab ' ITPC'], 'interpreter', 'none');
    set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',ylims);
    xlabel('Time (s)');
    ylabel('ITPC');
    legend([t_line.mainLine b_line.mainLine],{...
        ['Theta'],...% (' num2str(mean(theta_lim(:,ch_ix,1))) '-' num2str(mean(theta_lim(:,ch_ix,2))) ' Hz)
        ['Beta']},'Location','best');% (' num2str(mean(beta_lim(:,ch_ix,1))) '-' num2str(mean(beta_lim(:,ch_ix,2))) ' Hz)
    set(gca,'FontSize',font_size);
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

