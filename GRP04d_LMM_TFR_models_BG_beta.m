%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40';%'TFRmth_D1t1_zS8t0_f2t40';%
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
outlier_thresh = 4;
n_quantiles = 5;

save_fig = 1;
fig_ftype = 'png';

SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};

sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
out_thresh_str = ['_out' num2str(outlier_thresh)];

table_name = [an_id norm_bhv_str norm_nrl_str];
fig_dir   = [prj_dir 'results/TFR/LMM/' table_name out_thresh_str '/BG_betalo/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

%% Load data
bhvs       = cell(size(SBJs));
mdls       = cell(size(SBJs));
for s = 1:length(SBJs)
    % Load behavior
    load([prj_dir 'data/' SBJs{s} '/' SBJs{s} '_stim_preproc.mat'],'sbj_data');
    bhvs{s} = sbj_data.bhv;
    mdls{s} = sbj_data.mdl;
end

%% Load group model table
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Toss outliers
pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(pow_vars{f}) = table_all(~out_idx_all.(pow_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',pow_vars{f});
    if any(out_idx_all.(pow_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(pow_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(pow_vars{f})),pow_vars{f});
        fprintf(2,'%.2f, ',table_all.(pow_vars{f})(out_idx_all.(pow_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        pow_vars{f},size(good_tbl_all.(pow_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create current and previous trial tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
for p = 1:length(pow_vars)
    prv_nan_idx = isnan(good_tbl_prv.(pow_vars{p}).SV_prv);
    good_tbl_prv.(pow_vars{p})(prv_nan_idx,:) = [];
    prv_fields = good_tbl_prv.(pow_vars{p}).Properties.VariableNames;
    for f = 1:length(prv_fields)
        if any(isnan(good_tbl_prv.(pow_vars{p}).(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
    end
end

%% ========================================================================
%   BASAL GANGLIA BETA LOW
%  ========================================================================
lme_all = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + BG_roi + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ SV_cur + SV_prv + BG_roi + (1|sbj_n)');

%% Plot BG beta low by ROI
bg_roi_idx_all = good_tbl_all.BG_betalo.BG_roi;
figure;
errorbar([1 2],[mean(good_tbl_all.BG_betalo.BG_betalo(bg_roi_idx_all==1)) ...
    mean(good_tbl_all.BG_betalo.BG_betalo(bg_roi_idx_all==2))], [std(good_tbl_all.BG_betalo.BG_betalo(bg_roi_idx_all==1)) ...
    std(good_tbl_all.BG_betalo.BG_betalo(bg_roi_idx_all==2))]);
xlabel('BG ROI');
ylabel('BG beta low');
xlim([0.5 2.5]);
set(gca,'FontSize',16);

% BG modeling section from _orig:
% BG beta low and subjective value:
lme0 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG beta low and effort:
lme0 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ effort_cur + BG_roi + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ effortS_cur + BG_roi + (1|sbj_n)');
bg_betalo_eff = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
bg_betalo_effS = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Plot BG betalo ~ effortS as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_all.BG_betalo,'effortS_cur','BG_betalo',lme2,bg_betalo_effS.pValue(2));
xlabel('Subjective Effort (z)');
ylabel('BG low beta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot BG betalo ~ effortS as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.BG_betalo,'effortS_cur','BG_betalo',...
    lme2,bg_betalo_effS.pValue(2),n_quantiles);
xlabel('Subjective Effort (z)');
ylabel('BG low beta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% BG previous trial models
% BG beta low and subjective value:
lme0 = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ SV_prv + BG_roi + (1|sbj_n)');
bg_betalo_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG beta low and effort:
lme0 = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ effort_prv + BG_roi + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ effortS_prv + BG_roi + (1|sbj_n)');
bg_betalo_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
bg_betalo_effS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Plot BG betalo ~ previous subjective effort as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.BG_betalo,'effortS_prv','BG_betalo',lme2,bg_betalo_effS_prv.pValue(2));
xlabel('Prev. Subjective Effort (z)');
ylabel('BG low beta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


