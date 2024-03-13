%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'; stat_id = 'S5t15_bhvz_nrl0_out3';%_rt21';
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

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/PFC_beta/'];
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
pow_vars = {'PFC_betalo'};
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
%   PFC BETA MODELS
%  ========================================================================
%% Full Reward-Effort model
lme_full = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewc = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewp = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)');
lme_full_noeffc = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_noeffp = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');

pfc_betalo_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
pfc_betalo_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
pfc_betalo_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
pfc_betalo_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

lme_full_nop = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + effortS_cur + (1|sbj_n)');
pfc_betalo_addprv = compare(lme_full_nop,lme_full,'CheckNesting',true)

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ SV_cur + SV_prv + (1|sbj_n)');
pfc_betalo_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme_sv_cur = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ SV_cur + (1|sbj_n)');
lme_sv_prv = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ SV_prv + (1|sbj_n)');
pfc_betalo_svc = compare(lme_sv_prv,lme_sv_curprv,'CheckNesting',true)
pfc_betalo_svp = compare(lme_sv_cur,lme_sv_curprv,'CheckNesting',true)

% Plot PFC beta ~ SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_betalo,'SV_cur','PFC_betalo',...
    lme_sv_curprv,pfc_betalo_svc.pValue(2),n_quantiles,'continuous',1);
xlabel('Subjective Value (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('PFC beta (z)');
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

%% PFC Beta and Reward vs. Effort models:
% PFC beta low and reward:
lme0 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ reward_cur + (1|sbj_n)');
pfc_betalo_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% PFC beta low and effort:
lme0 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ effort_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ effortS_cur + (1|sbj_n)');
lme3 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ SV_cur + (1|sbj_n)');
pfc_betalo_eff = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_effS = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_SVc = compare(lme0,lme3,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_eff_vs_effS = compare(lme1,lme2,'NSim',1000)
pfc_betalo_effS_vs_SV = compare(lme3,lme2,'NSim',1000)
%   effortS is better model than effort
%   effort and effort S are better models than SV

% Plot PFC betalo ~ effortS as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_all.PFC_betalo,'effortS_cur','PFC_betalo',lme2,pfc_betalo_effS.pValue(2));
xlabel('Subjective Effort (z)');
ylabel('PFC low beta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC betalo ~ effortS as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.PFC_betalo,'effortS_cur','PFC_betalo',...
    lme2,pfc_betalo_effS.pValue(2),n_quantiles);
xlabel('Subjective Effort (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('PFC low beta (z)');
ylim([-0.6 0.6]);
yticks(-0.5:0.5:0.5);
set(gca,'FontSize',20);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC betalo ~ effortS as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_betalo,'effortS','PFC_betalo');
ylabel('PFC low beta (z)');
set(gca,'XTickLabel',{'Low Previous EFF','High Previous EFF'});
legend({'Low Current EFF','High Current EFF'},'Location','best');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% PFC beta low and previous effort:
lme0 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ effort_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ effortS_prv + (1|sbj_n)');
pfc_betalo_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_effS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

%% PFC beta and decision:
lme0 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ decision_cur + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.PFC_betalo,'PFC_betalo~ decision_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betalo_dec_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_dec_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

lme0 = fitglme(good_tbl_prv.PFC_betalo,'decision_cur ~ 1 + (1|sbj_n)','Distribution','binomial');
lme1 = fitglme(good_tbl_prv.PFC_betalo,'decision_cur ~ PFC_betalo + (1|sbj_n)','Distribution','binomial');
dec_cur_pfc_betalo = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% PFC beta and reward change and Global Reward State (GRS):
lme0 = fitlme(good_tbl_grs.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_grs.PFC_betalo,'PFC_betalo~ reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_grs.PFC_betalo,'PFC_betalo~ grs + (1|sbj_n)');
lme3 = fitlme(good_tbl_grs.PFC_betalo,'PFC_betalo~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betalo_rew_chg = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_grs     = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)


