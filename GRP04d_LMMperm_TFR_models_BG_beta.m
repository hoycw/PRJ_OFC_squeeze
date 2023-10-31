%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4_bt1k';%'S5t15_bhvz_nrl0_out4';%
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4';
% Pre-decision:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'Dn5t0_bhvz_nrlz_out4';
% Post-decision/feedback:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'D0t1_bhvz_nrlz_out4';% stat_id = 'D0t5_bhvz_nrlz_out4';%

n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if ~isfield(st,'nboots'); error('only run this for permutation bootstrap stat_ids!'); end

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/BG_beta/'];
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
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Toss outliers
pow_vars = {'BG_betalo'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>mean(table_all.(pow_vars{f}))...
                                    +(st.outlier_thresh*std(table_all.(pow_vars{f})));
    
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
        fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},st.outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        pow_vars{f},size(good_tbl_all.(pow_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create previous trial and GRS tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
tbl_fields = good_tbl_all.(pow_vars{1}).Properties.VariableNames;
for p = 1:length(pow_vars)
    prv_nan_idx = isnan(good_tbl_prv.(pow_vars{p}).SV_prv);
    good_tbl_prv.(pow_vars{p})(prv_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_prv.(pow_vars{p}).(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
            error(['NaN is table_prv.' tbl_fields{f}]);
        end
    end
end

% Toss NaNs from grs table
good_tbl_grs = good_tbl_all;
for p = 1:length(pow_vars)
    grs_nan_idx = isnan(good_tbl_grs.(pow_vars{p}).grs);
    good_tbl_grs.(pow_vars{p})(grs_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_grs.(pow_vars{p}).(tbl_fields{f})))
            fprintf(['%d NaNs in good_tbl_grs.' tbl_fields{f}]);
        end
    end
end

%% ========================================================================
%   BASAL GANGLIA BETA LOW
%  ========================================================================
%% Full reward/effort model
fit_method = 'ML';
full_formula = 'BG_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)');
predictors = {'reward_cur','effortS_cur','reward_prv','effortS_prv'};
[lme_full,pred_pval] = fn_run_LMM_full_permutation_models(good_tbl_prv.BG_betalo,predictors,full_formula,fit_method,st.nboots);

%% Plot BG beta ~ current effort as line plot
fit_method = 'ML';
full_formula = 'BG_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';
lme_full     = fitlme(good_tbl_prv.BG_betalo,full_formula,'FitMethod',fit_method);
noEc_formula = 'BG_betalo~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)';
lme_noRc     = fitlme(good_tbl_prv.BG_betalo,noEc_formula,'FitMethod',fit_method);
noEc_p = compare(lme_noRc,lme_full,'CheckNesting',true);

fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.BG_betalo,'effortS_cur','BG_betalo',...
    lme_full,noEc_p.pValue(2),n_quantiles);
xlabel('Subjective Effort (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('BG beta (z)');
if strcmp(st.norm_nrl_pred,'zscore')
    ylim([-0.6 0.6]);
    yticks(-0.5:0.5:0.5);
end
set(gca,'FontSize',20);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compare reward + effort vs. SV
fit_method = 'ML';
lme_all = fitlme(good_tbl_prv.BG_betalo,'BG_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');

sv_formula = 'BG_betalo~ SV_cur + SV_prv + (1|sbj_n)';
sv_pred = {'SV_cur','SV_prv'};
[lme_sv,sv_pred_pval] = fn_run_LMM_full_permutation_models(good_tbl_prv.BG_betalo,sv_pred,sv_formula,fit_method,st.nboots);
bg_beta_full_vs_SV = compare(lme_sv,lme_all,'NSim',1000)
