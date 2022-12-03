%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_zS8t0_f2t40';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'logz';%'none';%
outlier_thresh = 4;

save_fig = 1;
fig_ftype = 'png';

SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

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
if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
% out_thresh_str = ['out' num2str(outlier_thresh)];
table_name = [an_id norm_bhv_str norm_nrl_str];

% combined trial table
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Create current and previous trial tables
% Toss NaNs from previous table
table_prv = table_all;
prv_nan_idx = isnan(table_prv.SV_prv);
table_prv(prv_nan_idx,:) = [];
prv_fields = table_prv.Properties.VariableNames;
for f = 1:length(prv_fields)
    if any(isnan(table_prv.(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
end

%% Toss outliers
pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
out_idx_prv = struct;
out_idx_all = struct;
out_ix_prv = [];
out_ix_all = [];
good_tbl_prv = struct;
good_tbl_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_prv.(pow_vars{f}) = abs(table_prv.(pow_vars{f}))>outlier_thresh;
    out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_prv.(pow_vars{f}) = table_prv(~out_idx_prv.(pow_vars{f}),:);
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
    
    if any(out_idx_prv.(pow_vars{f}))
        out_ix_prv = [out_ix_prv; find(out_idx_prv.(pow_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_prv:\t',sum(out_idx_prv.(pow_vars{f})),pow_vars{f});
        fprintf(2,'%.2f, ',table_prv.(pow_vars{f})(out_idx_prv.(pow_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_prv = %d / %d\n',...
        pow_vars{f},size(good_tbl_prv.(pow_vars{f}),1),size(table_prv,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));
all_outliers_prv = unique(out_ix_prv);
fprintf(2,'Total bad trials in table_prv: %d\n',length(all_outliers_prv));

%% ========================================================================
%   BASAL GANGLIA MODELS
%  ========================================================================

% BG beta low and subjective value:
lme0 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ reward_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% BG modeling section from _orig:
fig_dir   = [prj_dir 'results/TFR/LMM/' table_name '/BG_betalo/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

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

% Plot betalo ~ current subjective effort as scatter plot
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_LFP_betalo_effortS_scatter';
figure('Name',fig_name); hold on;
xvals = min(good_tbl_all.BG_betalo.effortS_cur)-x_fudge:0.01:max(good_tbl_all.BG_betalo.effortS_cur)+x_fudge;
for s = 1:length(SBJs)
    scatter(good_tbl_all.BG_betalo.effortS_cur(good_tbl_all.BG_betalo.sbj_n==s),good_tbl_all.BG_betalo.BG_betalo(good_tbl_all.BG_betalo.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(good_tbl_all.BG_betalo.effortS_cur(good_tbl_all.BG_betalo.sbj_n==s),good_tbl_all.BG_betalo.BG_betalo(good_tbl_all.BG_betalo.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel('Subj. Effort (z)');
ylabel('Basal Ganglia low beta (z)');
title(['LMM coefficient = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(bg_betalo_effS.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
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

% Plot betalo ~ previous subjective effort as scatter plot
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_LFP_betalo_effortS_scatter';
figure('Name',fig_name); hold on;
xvals = min(good_tbl_prv.BG_betalo.effortS_prv)-x_fudge:0.01:max(good_tbl_prv.BG_betalo.effortS_prv)+x_fudge;
for s = 1:length(SBJs)
    scatter(good_tbl_prv.BG_betalo.effortS_prv(good_tbl_prv.BG_betalo.sbj_n==s),good_tbl_prv.BG_betalo.BG_betalo(good_tbl_prv.BG_betalo.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(good_tbl_prv.BG_betalo.effortS_prv(good_tbl_prv.BG_betalo.sbj_n==s),good_tbl_prv.BG_betalo.BG_betalo(good_tbl_prv.BG_betalo.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel('Subj. Effort (z)');
ylabel('Basal Ganglia low beta (z)');
title(['LMM coefficient = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(bg_betalo_effS_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% BG theta and previous subjective value:
lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + BG_roi + (1|sbj_n)');
bg_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
% lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + (1|sbj_n)');
% bg_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

