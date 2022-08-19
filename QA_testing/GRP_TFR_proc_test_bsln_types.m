%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
% clear all
% close all
% clc

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_ids     = {'TFRw_S25t201_z25t05_fl2t40_varWid'};%'simon'};%'TFRw_S25t201_z25t05_fl2t40','TFRm_S25t201_zbtS_sm0_l0_wnVar','simon'};
bsln_types = {'relchange','db','zscore','zboot'};
bsln_boots = 500;

save_fig  = 1;
plt_id    = 'ts_S2t2_evnts_sigLine';
fig_ftype = 'png';
fig_vis   = 'on';
fig_dir   = [prj_dir 'results/TFR/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

outlier_std_thresh = 3;

plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Time Frequency analysis
for s = 1:4
    %% Load data
    preproc_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [preproc_dir SBJs{s} '_full_preproc.mat'];
    fprintf('Loading %s\n',proc_fname);
    load(proc_fname,'sbj_data');
    
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
            
            % Re-align to decision (if needed)
%             [tfr_r] = fn_realign_tfr_s2r(tfr_s,rt,trial_lim_sec);
            
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
            'outerposition',[0 0 1 1],'Visible',fig_vis);
        
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
        
        % Save Figure
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
        
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
    
end

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
