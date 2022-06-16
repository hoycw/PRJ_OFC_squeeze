%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
clc

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_ids = {'TFRw_S25t201_z25t05_fl2t40','TFRm_S25t201_zbtS_sm0_l0_wnVar'};
bsln_types = {'relchange','zscore','zboot'};
bsln_boots = 500;
plt_id     = 'ts_S2t2_evnts_sigLine';

outlier_std_thresh = 3;

plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Time Frequency analysis
for s = 1:4
    %% Load data
    preproc_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [preproc_dir SBJs{s} '_' an.evnt_lab '_preproc.mat'];
    fprintf('Loading %s\n',proc_fname);
    data{s} = load(proc_fname,'sbj_data');
    
    %% Time-frequency representation
    tfr = cell([numel(an_ids) numel(bsln_types)]);
    for an_ix = 1:length(an_ids)
        an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_ids{an_ix} '_vars.m'];
        eval(an_vars_cmd);
        
        % Trim to analysis time
        [trial_lim_s_pad] = fn_get_filter_padding(cfg_tfr,an.trial_lim_s,an.bsln_lim);
        cfgs = []; cfgs.latency = trial_lim_s_pad;
        trim = ft_selectdata(cfgs,sbj_data.ts);
        
        % Extract time-frequency representation
        tfr_raw = ft_freqanalysis(cfg_tfr, trim);
        
        % Baseline correction
        for b_ix = 1:numel(bsln_types)
            switch bsln_types{b_ix}
                case {'zscore', 'zboot', 'demean', 'my_relchange'}
                    tfr_bsln = fn_bsln_ft_tfr(tfr_raw,an.bsln_lim,bsln_types{b_ix},bsln_boots);
                case {'relchange','db'}
                    cfgbsln = [];
                    cfgbsln.baseline     = an.bsln_lim;
                    cfgbsln.baselinetype = bsln_types{b_ix};
                    cfgbsln.parameter    = 'powspctrm';
                    tfr_bsln = ft_freqbaseline(cfgbsln,tfr_raw);
                case 'none'
                    if an.complex
                        fprintf('\tNo baseline correction for ITPC data...\n');
                    else
                        error('Why no baseline correction if not ITPC?');
                    end
                otherwise
                    error(['No baseline implemented for bsln_type: ' bsln_types{b_ix}]);
            end
            
            % Trim back down to analysis window
            cfgs = [];
            % cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
            cfgs.channel = an.ROI;
            cfgs.latency = an.trial_lim_s;
            cfgs.avgoverrpt = 'yes';
            tfr{an_ix,b_ix} = ft_selectdata(cfgs, tfr_bsln);
        end
    end
    
    
    %% Plot TFRs
    for an_ix = 1:numel(an_ids)
        fig_name = [SBJs{s} '_' an_ids{an_ix}];
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
        
        for ch_ix = 1:numel(tfr{an_ix,b_ix}.label)
            for b_ix = 1:numel(bsln_types)
                subplot(length(tfr{an_ix,b_ix}.label),length(bsln_types),ch_ix*length(bsln_types)-b_ix+1);
                % Get color lims per condition
                vals = tfr{an_ix,b_ix}.powspctrm(ch_ix,:,:);
                clim = [min(vals(:)) max(vals(:))];
                
                % Plot TFR
                imagesc(tfr{an_ix,b_ix}.time, tfr{an_ix,b_ix}.freq, squeeze(tfr{an_ix,b_ix}.powspctrm(ch_ix,:,:)));
                set(gca,'YDir','normal');
                caxis(clim);
                colorbar;
                
                % Plot Events
                line([0 0],ylim,...
                    'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                    'LineStyle',plt.evnt_styles{1});
                
                % Axes and parameters
                title([tfr{an_ix,b_ix}.label{ch_ix} '- ' bsln_types{b_ix}]);
                set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
                set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                set(gca,'FontSize',16);
            end
        end
    end
    
    %%
    if plot_it
        %     cfg = [];
        %     figure;
        %     subplot(2,1,1);
        %     cfg.channel      = 'OFC';
        %     ft_singleplotTFR(cfg,freqout);
        %     xlim([-1.5 1.5]);
        %     subplot(2,1,2);
        %     cfg.channel      = 'LFP';
        %     ft_singleplotTFR(cfg,freqout);
        %     xlim([-1.5 1.5]);
        
        cfg = [];
        figure;
        subplot(2,1,1);
        cfg.channel      = 'OFC';
        ft_singleplotTFR(cfg,sbj_data.TFbl);
        xlim([-0.5 2]);
        ylim([ 2 35]);
        subplot(2,1,2);
        cfg.channel      = 'LFP';
        ft_singleplotTFR(cfg,sbj_data.TFbl);
        xlim([-0.5 2]);
        ylim([ 2 35]);
    end
    
    %% SAVE data
%     sbj_data.desp    = preproc_name;
%     sbj_data.ts      = data;
%     sbj_data.TF      = freqout;
%     sbj_data.TFbl    = freqoutbl;
%     sbj_data.bhv     = bhv;
%     sbj_data.mdl.par = par;
%     sbj_data.mdl.fit = fit;
    
%     preproc_dir = [prj_dir 'data/' SBJs{s} '/'];
%     proc_fname = [preproc_dir SBJs{s} '_' evnt_lab '_preproc.mat'];
%     fprintf('Loading %s\n',proc_fname);
%     load(proc_fname,'sbj_data');
    
    %% Statistical testing %%
    
    % Start by examining for a change in the power with the event %
    
    % ft_statistics_montecarlo
    
    % Start by examining the reward and effort separately %
    % compute statistics with ft_statfun_indepsamplesregrT
    
%     cfg = [];
%     cfg.channel          = 'OFC';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     % cfg.statistic        = 'ft_statfun_depsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     % This should be 1 - 3 as this actually include the post presentation
%     % period only.
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = stake';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statOFCstake = ft_freqstatistics(cfg, freqoutbl);
%     
%     cfg = [];
%     cfg.channel          = 'OFC';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     % cfg.statistic        = 'ft_statfun_depsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = EFF(par(2))';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statOFCeffort = ft_freqstatistics(cfg, freqoutbl);
%     
%     cfg = [];
%     cfg.channel          = 'OFC';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     % cfg.statistic        = 'ft_statfun_depsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = SubjVal1';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statOFCSV = ft_freqstatistics(cfg, freqoutbl);
%     
%     cfg = [];
%     cfg.channel          = 'LFP';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = stake';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statLFPstake = ft_freqstatistics(cfg, freqoutbl);
%     
%     cfg = [];
%     cfg.channel          = 'LFP';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = EFF(par(2))';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statLFPeffort = ft_freqstatistics(cfg, freqoutbl);
%     
%     cfg = [];
%     cfg.channel          = 'LFP';
%     cfg.statistic        = 'ft_statfun_indepsamplesregrT';
%     cfg.method           = 'montecarlo';
%     cfg.correctm         = 'cluster';
%     cfg.numrandomization = 1000;
%     cfg.alpha            = 0.05;
%     cfg.tail             = 0;
%     cfg.correcttail      ='alpha';
%     cfg.frequency        = [2 30];
%     cfg.latency          = [0 1.5];
%     
%     n1 = size(stake,1);
%     design(1,1:n1)       = SubjVal1';
%     
%     cfg.design           = design;
%     cfg.ivar             = 1;
%     
%     statLFPSV = ft_freqstatistics(cfg, freqoutbl);
%     
%     % PLOT STATS RESULTS
%     
%     figure;
%     subplot(2,3,1);
%     cfg           = [];
%     cfg.channel   = {'OFC'};
%     cfg.parameter = 'stat';
%     cfg.colormap  = parula;
%     cfg.ylim      = [2 100];
%     cfg.xlim      = [-1.5 1.5];
%     % cfg.zlim      = [ ] ;
%     cfg.marker    ='off';
%     cfg.style     = 'fill';
%     cfg.comment   = 'off';
%     cfg.maskparameter = 'mask';
%     cfg.maskstyle = 'outline';
%     cfg.colorbar  = 'yes';
%     ft_singleplotTFR(cfg,statOFCstake);
%     
%     subplot(2,3,2);
%     cfg.channel   = {'OFC'};
%     ft_singleplotTFR(cfg,statOFCeffort);
%     
%     subplot(2,3,3);
%     cfg.channel   = {'OFC'};
%     ft_singleplotTFR(cfg,statOFCSV);
%     
%     subplot(2,3,4);
%     cfg.channel   = {'LFP'};
%     ft_singleplotTFR(cfg,statLFPstake);
%     
%     subplot(2,3,5);
%     cfg.channel   = {'LFP'};
%     ft_singleplotTFR(cfg,statLFPeffort);
%     
%     subplot(2,3,6);
%     cfg.channel   = {'LFP'};
%     ft_singleplotTFR(cfg,statLFPSV);
end