%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
% sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_ids = {'TFRmth_S1t2_dbS1t0_f2t40','TFRmth_S1t2_zS1t0_f2t40','TFRmth_S1t2_zS25t05_f2t40'};%'TFRw_S25t2_dbS25t05_fl2t40_c7','TFRw_D1t1_dbS25t05_fl2t40_c7'};
%'TFRw_S25t2_noBsln_fl1t40_c7','TFRw_S25t2_zbtS25t05_fl1t40_c7'};%'TFRw_S25t2_noBsln_fl2t40_c7'};%

if contains(an_ids{1},'_S')
    an_lim = [0.5 1.5];
elseif contains(an_ids{1},'_D')
    an_lim = [-0.5 0];
end

freq_names = {'theta','beta low'};
freq_bands = [4 7; 12 20];
mrkr_sz = 25;

%!!! add checks for same TFR method and same time lock event

%% Time Frequency analysis
an_vars = cell(size(an_ids));
tfr_raw = cell([numel(SBJs) numel(an_ids)]);
tfr_trim = cell([numel(SBJs) numel(an_ids)]);
tfr_bsln = cell([numel(SBJs) numel(an_ids)]);
tfr = cell([numel(SBJs) numel(an_ids)]);
for s = 1:4
    %% Load data
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    tfr_fname = [sbj_dir SBJs{s} '_stim_preproc.mat']; % This is [-3 10] to cover all possible times
    fprintf('Loading %s\n',tfr_fname);
    load(tfr_fname,'sbj_data');
    
    %% Time-frequency representation
    for an_ix = 1:length(an_ids)
        an_vars_cmd = ['run ' prj_dir 'scripts/an_vars/' an_ids{an_ix} '_vars.m'];
        eval(an_vars_cmd);
        an_vars{an_ix} = an;
        
%         % Trim to analysis time
%         if strcmp(an.event_type,'S')
%             [trial_lim_s_pad] = fn_get_filter_padding(cfg_tfr,an.trial_lim_s,an.bsln_lim);
%             cfgs = []; cfgs.latency = trial_lim_s_pad;
%             trim = ft_selectdata(cfgs,sbj_data.ts);
%         else
%             % That function doesn't work with variable length trials
%             trim = sbj_data.ts;
%         end
        
        % Extract time-frequency representation
        tfr_raw{s,an_ix} = ft_freqanalysis(cfg_tfr, sbj_data.ts);%trim);
        
        % Trim back down to original trial_lim_s to exclude NaNs
        cfg_trim = [];
        if strcmp(an.event_type,'S')
            cfg_trim.latency = an.trial_lim_s;
        elseif strcmp(an.event_type,'D') && strcmp(an.bsln_evnt,'S')
            cfg_trim.latency = [an.bsln_lim(1) max(sbj_data.bhv.rt)+an.trial_lim_s(2)];
        else
            error('mismatched event without S-locked baseline!');
        end
        tfr_trim{s,an_ix} = ft_selectdata(cfg_trim,tfr_raw{s,an_ix});
        
        % Extract raw baseline
        cfg_bsln = [];
        cfg_bsln.latency = an.bsln_lim;
        tfr_bsln{s,an_ix} = ft_selectdata(cfg_bsln,tfr_raw{s,an_ix});
        
        % Baseline correction
        switch an.bsln_type
            case {'zscore','zboot', 'demean', 'my_relchange'}
                tfr{s,an_ix} = fn_bsln_ft_tfr(tfr_trim{s,an_ix},an.bsln_lim,an.bsln_type,an.bsln_boots);
            case {'relchange','db'}
                cfgbsln = [];
                cfgbsln.baseline     = an.bsln_lim;
                cfgbsln.baselinetype = an.bsln_type;
                cfgbsln.parameter    = 'powspctrm';
                tfr{s,an_ix} = ft_freqbaseline(cfgbsln,tfr_trim{s,an_ix});
            case 'none'
                warning('\tSkipping baseline correction!');
                tfr{s,an_ix} = tfr_trim{s,an_ix};
            otherwise
                error(['No baseline implemented for bsln_type: ' an.bsln_type]);
        end
        
        % Re-align to decision and re-trim
        if ~strcmp(an.event_type,an.bsln_evnt)
            if strcmp(an.event_type,'D') && strcmp(an.bsln_evnt,'S')
                [tfr{s,an_ix}] = fn_realign_tfr_s2r(tfr{s,an_ix},sbj_data.bhv.rt,an.trial_lim_s);
            else
                error('unknown combination of analysis and baseline events');
            end
        end
    end
end

%% Plot baseline values over time
an_ids_str = strjoin(an_ids,'-');
% [n_rc,~] = fn_num_subplots(numel(an_ids));

ch_ix = 2;
for s = 1:4
    trl_colors = jet(size(tfr_bsln{s,1}.powspctrm,1));
    fig_name = [SBJs{s} '_TFR_bsln_trial_scatter_' an_ids_str];
    figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.6 0.5]);
    for band_ix = 1:size(freq_bands,1)
        for an_ix = 1:numel(an_ids)
            % Extract power values
            cfg = [];
            cfg.avgoverrpt = 'no';
            cfg.avgoverfreq = 'yes';
            cfg.avgovertime = 'yes';
            cfg.frequency = freq_bands(band_ix,:);
            bsln = ft_selectdata(cfg,tfr_bsln{s,an_ix});
            
            % Plot over time
            subplot(size(freq_bands,1),numel(an_ids),band_ix*numel(an_ids)-an_ix+1);
            scatter(1:size(bsln.powspctrm,1),bsln.powspctrm(:,ch_ix),mrkr_sz,trl_colors);
            xlabel('Trial');
            ylabel('Baseline Power');
            title([freq_names{band_ix} ': ' an_ids{an_ix}],'Interpreter','none');
            set(gca,'FontSize',16);
        end
    end
end

%% Plot baseline vs. main effect
ch_ix = 2;
for s = 2:4
    trl_colors = jet(size(tfr_raw{s,1}.powspctrm,1));
    fig_name = [SBJs{s} '_TFR_bsln_vs_trial_violins_' an_ids_str];
    figure('Name',fig_name,'units','norm','OuterPosition',[0 0 1 1]);
    for band_ix = 1:size(freq_bands,1)
        for an_ix = 1:numel(an_ids)
            % Extract power values
            cfg = [];
            cfg.avgoverrpt  = 'no';
            cfg.avgoverfreq = 'yes';
            cfg.avgovertime = 'yes';
            cfg.frequency   = freq_bands(band_ix,:);
            cfg.latency     = an_vars{an_ix}.bsln_lim;
            bsln = ft_selectdata(cfg,tfr_raw{s,an_ix});
            
            cfg.latency     = an_lim;
            trl = ft_selectdata(cfg,tfr_raw{s,an_ix});
            
            % Plot distribution of values
            subplot(size(freq_bands,1),numel(an_ids),band_ix*numel(an_ids)-an_ix+1);
            violins = violinplot([bsln.powspctrm(:,ch_ix) trl.powspctrm(:,ch_ix)],{'Baseline','Trial'},'ShowMean',true);%,'ViolinAlpha',0.3);
            
            % Adjust plot propeties
            for violin_ix = 1:2
                % Fix mean line
                violins(violin_ix).MeanPlot.LineWidth = 3;
%                 % Violin mean line is off since ksdensity is tossing some of the
%                 % values for some reason, so manually set it and adjust width
%                 violins(violin_ix).MeanPlot.YData = repmat(mean(plot_onsets{reg_ix}.(roi_list{roi_ix})),1,2);
%                 mean_ix = nearest(violins(violin_ix).ViolinPlot.YData,mean(plot_onsets{reg_ix}.(roi_list{roi_ix})));
%                 violins(violin_ix).MeanPlot.XData = [violin_ix-(violins(violin_ix).ViolinPlot.XData(mean_ix)-violin_ix), ...
%                     violins(violin_ix).ViolinPlot.XData(mean_ix)];
%                 violins(violin_ix).ViolinColor = [0.8 0.8 0.8];
%                 violins(violin_ix).BoxPlot.FaceColor = roi_colors{roi_ix};
%                 violins(violin_ix).EdgeColor = roi_colors{roi_ix};
                
                % Change scatter colors within violin to mark trial
                violins(violin_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
                violins(violin_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
                violins(violin_ix).ScatterPlot.CData = trl_colors;
            end
            
            % Draw lines between same trial
            decrease_idx = violins(1).ScatterPlot.YData > violins(2).ScatterPlot.YData;
            for t_ix = 1:size(tfr_raw{s,an_ix}.powspctrm,1)
                if decrease_idx(t_ix)
                    line([violins(1).ScatterPlot.XData(t_ix) violins(2).ScatterPlot.XData(t_ix)],...
                        [violins(1).ScatterPlot.YData(t_ix) violins(2).ScatterPlot.YData(t_ix)],...
                        'Color','b');
                else
                    line([violins(1).ScatterPlot.XData(t_ix) violins(2).ScatterPlot.XData(t_ix)],...
                        [violins(1).ScatterPlot.YData(t_ix) violins(2).ScatterPlot.YData(t_ix)],...
                        'Color','r');
                end
            end
            
            ylabel('Power');
            title([freq_names{band_ix} ': ' an_ids{an_ix}],'Interpreter','none');
            set(gca,'FontSize',16);
        end
    end
end

%% Compare proportion of increases and decreases between raw and baseline-corrected trials
% histogram of trl - bsln for raw and all an_ids
n_bins = 50;
for s = 1:4
    trl_colors = jet(size(tfr_raw{s,1}.powspctrm,1));
    fig_name = [SBJs{s} '_TFR_bsln_trl_change_hist_' an_ids_str];
    figure('Name',fig_name,'units','norm','OuterPosition',[0 0 1 1]);
    for band_ix = 1:size(freq_bands,1)
        for an_ix = 1:numel(an_ids)+1
            % Extract power values
            cfg = [];
            cfg.avgoverrpt  = 'no';
            cfg.avgoverfreq = 'yes';
            cfg.avgovertime = 'yes';
            cfg.frequency   = freq_bands(band_ix,:);
            if an_ix==1
                cfg.latency     = an_vars{1}.bsln_lim;
                bsln = ft_selectdata(cfg,tfr_raw{s,1});
                cfg.latency     = an_lim;
                trl = ft_selectdata(cfg,tfr_raw{s,1});
            else
                cfg.latency     = an_vars{an_ix-1}.bsln_lim;
                bsln = ft_selectdata(cfg,tfr{s,an_ix-1});
                cfg.latency     = an_lim;
                trl = ft_selectdata(cfg,tfr{s,an_ix-1});
            end
            
            % Plot histogram of differences
            subplot(size(freq_bands,1),numel(an_ids)+1,band_ix*(numel(an_ids)+1)-an_ix+1);
            hold on;
            change = trl.powspctrm(:,ch_ix)-bsln.powspctrm(:,ch_ix);
            ch_hist   = histogram(change,n_bins,'FaceColor','r');
            trl_hist  = histogram(trl.powspctrm(:,ch_ix),n_bins,'FaceColor','g');
            bsln_hist = histogram(bsln.powspctrm(:,ch_ix),n_bins,'FaceColor','b');
            line([mean(change) mean(change)],ylim,'Color','r','LineWidth',2);
            
            xlabel('Power');%'Trial - Baseline Power');
%             xlim([-max(abs(change)) max(abs(change))]);
            legend([bsln_hist trl_hist ch_hist],{'Baseline','Trial','Change'},'Location','best');
            if an_ix==1
                title([freq_names{band_ix} ': Raw'],'Interpreter','none');
            else
                title([freq_names{band_ix} ': ' an_ids{an_ix-1}],'Interpreter','none');
            end
            set(gca,'FontSize',16);
        end
    end
end