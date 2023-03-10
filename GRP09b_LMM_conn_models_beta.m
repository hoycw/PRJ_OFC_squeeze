%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'ampcorr_Sn8t0_bhvz_nrlfz_out4';
% Stimulus decision phase:
% an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'ampcorr_S5t15_bhvz_nrlfz_out4';
% Pre-decision:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'ampcorr_Dn5t0_bhvz_nrlfz_out4';
% Post-decision/feedback:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'ampcorr_D0t5_bhvz_nrlfz_out4';

% Phase-locking value
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'PLV_S5t15_bhvz_nrlz_out4';

% Jackknife coherence
% an_id = 'TFRmth_S03t2_f2t30_fourier'; stat_id = 'cohjk_S5t15_bhvz_nrlz';

n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/beta_conn/'];
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
conn_vars = {'betalo_conn'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(conn_vars)
    % Identify outliers
    out_idx_all.(conn_vars{f}) = abs(table_all.(conn_vars{f}))>st.outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(conn_vars{f}) = table_all(~out_idx_all.(conn_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',conn_vars{f});
    if any(out_idx_all.(conn_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(conn_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(conn_vars{f})),conn_vars{f});
        fprintf(2,'%.2f, ',table_all.(conn_vars{f})(out_idx_all.(conn_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',conn_vars{f},st.outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        conn_vars{f},size(good_tbl_all.(conn_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create previous trial and GRS tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
tbl_fields = good_tbl_all.(conn_vars{1}).Properties.VariableNames;
for p = 1:length(conn_vars)
    prv_nan_idx = isnan(good_tbl_prv.(conn_vars{p}).SV_prv);
    good_tbl_prv.(conn_vars{p})(prv_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_prv.(conn_vars{p}).(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
            error(['NaN is table_prv.' tbl_fields{f}]);
        end
    end
end

% Toss NaNs from grs table
good_tbl_grs = good_tbl_all;
for p = 1:length(conn_vars)
    grs_nan_idx = isnan(good_tbl_grs.(conn_vars{p}).grs);
    good_tbl_grs.(conn_vars{p})(grs_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_grs.(conn_vars{p}).(tbl_fields{f})))
            fprintf(['%d NaNs in good_tbl_grs.' tbl_fields{f}]);
        end
    end
end

%% ========================================================================
%   LOW BETA CONNECTIVITY MODELS
%  ========================================================================
%% Full model
% lme_full = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEc = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEp = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewc = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewp = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + dec_ease_cur + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffc = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffp = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');
% 
% bconn_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% bconn_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% bconn_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true)
% bconn_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
% bconn_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
% bconn_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewc = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewp = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)');
lme_full_noeffc = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_noeffp = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');
lme_full_bg_roi = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + BG_roi + (1|sbj_n)');

bconn_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
bconn_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
bconn_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
bconn_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

bconn_full_bg_roi = compare(lme_full,lme_full_bg_roi,'CheckNesting',true)

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ SV_cur + SV_prv + (1|sbj_n)');
betalo_conn_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme_ez_curprv = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ dec_ease_cur + dec_ease_prv + (1|sbj_n)');

%% Check ROI effects
lme_bg_roi = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ BG_roi + (1|sbj_n)');
lme_pfc_roi = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ PFC_roi + (1|sbj_n)');

%% Low Beta connectivity and Reward vs. Effort models:
% Low beta connectivity low and reward:
lme0 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ reward_cur + (1|sbj_n)');
betalo_conn_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Low beta connectivity low and effort:
lme0 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ 1 + (1|sbj_n)');
% lme0bg = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ BG_roi + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ effort_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ effortS_cur + (1|sbj_n)');
lme3 = fitlme(good_tbl_all.betalo_conn,'betalo_conn~ SV_cur + (1|sbj_n)');
% betalo_conn_bgroi = compare(lme0,lme0bg,'CheckNesting',true)%,'NSim',1000)
betalo_conn_eff = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
betalo_conn_effS = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% betalo_conn_effS_vs_SV = compare(lme3,lme2,'NSim',1000)
%   effortS is better model than effort
%   effort and effort S are better models than SV

% Plot Low beta connectivity ~ effortS as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_all.betalo_conn,'effortS_cur','betalo_conn',lme2,betalo_conn_effS.pValue(2));
xlabel('Subjective Effort (z)');
ylabel('Low beta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot Low beta connectivity ~ effortS as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.betalo_conn,'effortS_cur','betalo_conn',...
    lme2,betalo_conn_effS.pValue(2),n_quantiles);
xlabel('Subjective Effort (z)');
ylabel('Low beta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot Low beta connectivity ~ effortS as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.betalo_conn,'effortS','betalo_conn');
ylabel('Low beta connectivity (z)');
set(gca,'XTickLabel',{'Low Previous EFF','High Previous EFF'});
legend({'Low Current EFF','High Current EFF'},'Location','best');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Low beta connectivity low and previous effort:
lme0 = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ effort_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.betalo_conn,'betalo_conn~ effortS_prv + (1|sbj_n)');
betalo_conn_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
betalo_conn_effS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

%% PFC beta lo and reward change and Global Reward State (GRS):
lme0 = fitlme(good_tbl_grs.betalo_conn,'betalo_conn~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_grs.betalo_conn,'betalo_conn~ reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_grs.betalo_conn,'betalo_conn~ grs + (1|sbj_n)');
lme3 = fitlme(good_tbl_grs.betalo_conn,'betalo_conn~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
betalo_conn_rew_chg = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
betalo_conn_grs     = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)


