%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc

addpath(['/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/']);
addpath(['/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/']);
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
% sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

% an_ids = {'TFRmth_S1t2_dbS1t0_f2t40','TFRmth_S1t2_zbtS1t0_f2t40','TFRmth_S1t2_zS1t0_f2t40','TFRmth_S1t2_zS25t05_f2t40'};%'TFRw_S25t2_dbS25t05_fl2t40_c7','TFRw_D1t1_dbS25t05_fl2t40_c7'};
%'TFRw_S25t2_noBsln_fl1t40_c7','TFRw_S25t2_zbtS25t05_fl1t40_c7'};%'TFRw_S25t2_noBsln_fl2t40_c7'};%
an_ids = {'TFRmth_S1t2_dbS1t0_f2t40','TFRmth_S1t2_zS1t0_f2t40','TFRmth_S1t2_zS1t0_f2t40_log','TFRmth_S1t2_zbtS1t0_f2t40'};%'TFRw_S25t2_dbS25t05_fl2t40_c7','TFRw_D1t1_dbS25t05_fl2t40_c7'};

if contains(an_ids{1},'_S')
    an_lim = [0.5 1.5];
elseif contains(an_ids{1},'_D')
    an_lim = [-0.5 0];
end

freq_names = {'theta','beta low'};
freq_bands = [4 7; 12 20];
ch_lab     = {'GPi','STN','GPi','STN'; 'FPC', 'OFC', 'OFC', 'FPC'};
mrkr_sz = 25;
save_fig = 1;
fig_ftype = 'png';

an_ids_str = strjoin(an_ids,'-');
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';

% Check for same TFR method and same time lock event
if ~all(contains(an_ids,'_S')) && ~all(contains(an_ids,'_D'))
    error('all TFRs should be locked to the same event!');
end
tfr_method = cell(size(an_ids));
for an_ix = 1:numel(an_ids)
    uscore = strfind(an_ids{an_ix},'_');
    tfr_method{an_ix} = an_ids{an_ix}(1:uscore(1)-1);
end
if ~all(strcmp(tfr_method{1},tfr_method))
    error('all TFRs should use same TFR method!');
end

%% Time Frequency analysis
an_vars = cell(size(an_ids));
tfr_raw = cell([numel(SBJs) numel(an_ids)]);
% tfr_bsln = cell([numel(SBJs) numel(an_ids)]);
tfr = cell([numel(SBJs) numel(an_ids)]);
for s = 1:4
    %% Load data
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    
    % Load Time-frequency representation
    for an_ix = 1:length(an_ids)
        an_vars_cmd = ['run ' prj_dir 'scripts/an_vars/' an_ids{an_ix} '_vars.m'];
        eval(an_vars_cmd);
        an_vars{an_ix} = an;
        
        tfr_fname = [sbj_dir SBJs{s} '_' an_ids{an_ix} '.mat'];
        fprintf('Loading %s\n',tfr_fname);
        tmp = load(tfr_fname);
        tfr_raw{s,an_ix}  = tmp.tfr_raw;
        tfr{s,an_ix} = tmp.tfr;
    end
end

%% Plot baseline values over time
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

%% Plot power distributions at baseline, trial, and change
% histogram of trl - bsln for raw and all an_ids
fig_dir   = [prj_dir 'results/TFR/bsln_properties/' an_ids_str '/bsln_trl_change/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
ch_ix = 2;
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
                title([ch_lab{ch_ix,s} ' ' freq_names{band_ix} ': Raw'],'Interpreter','none');
            else
                title([an_ids{an_ix-1}],'Interpreter','none');
            end
            set(gca,'FontSize',16);
        end
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot power distributions in baseline period
% histogram of trl - bsln for raw and all an_ids
fig_dir   = [prj_dir 'results/TFR/bsln_properties/' an_ids_str '/bsln/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
ch_ix = 2;
n_bins = 50;
for s = 1:4
    trl_colors = jet(size(tfr_raw{s,1}.powspctrm,1));
    fig_name = [SBJs{s} '_TFR_bsln_hist_' an_ids_str];
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
            else
                cfg.latency     = an_vars{an_ix-1}.bsln_lim;
                bsln = ft_selectdata(cfg,tfr{s,an_ix-1});
            end
            
            % Plot histogram of differences
            subplot(size(freq_bands,1),numel(an_ids)+1,band_ix*(numel(an_ids)+1)-an_ix+1);
            hold on;
            bsln_hist = histogram(bsln.powspctrm(:,ch_ix),n_bins,'FaceColor','b');
            line([mean(bsln.powspctrm(:,ch_ix)) mean(bsln.powspctrm(:,ch_ix))],ylim,'Color','r','LineWidth',2);
            
            xlabel('Baseline Power');%'Trial - Baseline Power');
%             xlim([-max(abs(change)) max(abs(change))]);
            if an_ix==1
                title([ch_lab{ch_ix,s} ' ' freq_names{band_ix} ': Raw'],'Interpreter','none');
            else
                title([an_ids{an_ix-1}],'Interpreter','none');
            end
            set(gca,'FontSize',16);
        end
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot power distributions in trial analysis period
% histogram of trl - bsln for raw and all an_ids
fig_dir   = [prj_dir 'results/TFR/bsln_properties/' an_ids_str '/trl/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
ch_ix = 2;
n_bins = 50;
for s = 1:4
    trl_colors = jet(size(tfr_raw{s,1}.powspctrm,1));
    fig_name = [SBJs{s} '_TFR_trl_hist_' an_ids_str];
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
                cfg.latency     = an_lim;
                trl = ft_selectdata(cfg,tfr_raw{s,1});
            else
                cfg.latency     = an_lim;
                trl = ft_selectdata(cfg,tfr{s,an_ix-1});
            end
            
            % Plot histogram of differences
            subplot(size(freq_bands,1),numel(an_ids)+1,band_ix*(numel(an_ids)+1)-an_ix+1);
            hold on;
            trl_hist  = histogram(trl.powspctrm(:,ch_ix),n_bins,'FaceColor','g');
            line([mean(trl.powspctrm(:,ch_ix)) mean(trl.powspctrm(:,ch_ix))],ylim,'Color','r','LineWidth',2);
            
            xlabel('Analysis Power');%'Trial - Baseline Power');
%             xlim([-max(abs(change)) max(abs(change))]);
            if an_ix==1
                title([ch_lab{ch_ix,s} ' ' freq_names{band_ix} ': Raw'],'Interpreter','none');
            else
                title([an_ids{an_ix-1}],'Interpreter','none');
            end
            set(gca,'FontSize',16);
        end
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end
