%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
close all
clear all

%%
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'; stat_id = 'S5t15_bhvz_nrl0_out3';

n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);

fig_dir   = [prj_dir 'results/bhv/LMM/' an_id '/' stat_id '/decision/'];
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

%% Load group model tables
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Toss outliers
dec_vars = {'decision_cur'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(dec_vars)
    % Identify outliers
    out_idx_all.(dec_vars{f}) = abs(table_all.(dec_vars{f}))>st.outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(dec_vars{f}) = table_all(~out_idx_all.(dec_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',dec_vars{f});
    if any(out_idx_all.(dec_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(dec_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(dec_vars{f})),dec_vars{f});
        fprintf(2,'%.2f, ',table_all.(dec_vars{f})(out_idx_all.(dec_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',dec_vars{f},st.outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        dec_vars{f},size(good_tbl_all.(dec_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create previous trial and GRS tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
tbl_fields = good_tbl_all.(dec_vars{1}).Properties.VariableNames;
for p = 1:length(dec_vars)
    prv_nan_idx = isnan(good_tbl_prv.(dec_vars{p}).SV_prv);
    good_tbl_prv.(dec_vars{p})(prv_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_prv.(dec_vars{p}).(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
            error(['NaN is table_prv.' tbl_fields{f}]);
        end
    end
end

% Toss NaNs from grs table
good_tbl_grs = good_tbl_all;
for p = 1:length(dec_vars)
    grs_nan_idx = isnan(good_tbl_grs.(dec_vars{p}).grs);
    good_tbl_grs.(dec_vars{p})(grs_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_grs.(dec_vars{p}).(tbl_fields{f})))
            fprintf(['%d NaNs in good_tbl_grs.' tbl_fields{f}]);
        end
    end
end

%% Create variables
% Mean Reward and Effort
efforts = unique(bhvs{1}.effort);
stakes  = unique(bhvs{1}.stake);
for s = 1:length(SBJs)
    if length(unique(bhvs{s}.effort))~=5 || length(unique(bhvs{s}.stake))~=5
        error([SBJs{s} ' is missing conditions for stake or effort!']);
    end
    bhvs{s}.SV_stake_mn  = nan([5 1]);
    bhvs{s}.SV_stake_se  = nan([5 1]);
    bhvs{s}.SV_effort_mn = nan([5 1]);
    bhvs{s}.SV_effort_se = nan([5 1]);
    bhvs{s}.dec_stake_mn  = nan([5 1]);
    bhvs{s}.dec_stake_se  = nan([5 1]);
%     bhvs{s}.dec_effort_mn = nan([5 1]);
%     bhvs{s}.dec_effort_se = nan([5 1]);
    for i = 1:5
        bhvs{s}.SV_effort_mn(i) = mean(bhvs{s}.SV(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_effort_se(i) = std(bhvs{s}.SV(bhvs{s}.effort==efforts(i)))./sqrt(sum(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_stake_mn(i) = mean(bhvs{s}.SV(bhvs{s}.stake==stakes(i)));
        bhvs{s}.SV_stake_se(i) = std(bhvs{s}.SV(bhvs{s}.stake==stakes(i)))./sqrt(sum(bhvs{s}.stake==stakes(i)));
%         bhvs{s}.dec_effort_mn(i) = mean(bhvs{s}.decision(bhvs{s}.effort==efforts(i)));
%         bhvs{s}.dec_effort_se(i) = std(bhvs{s}.decision(bhvs{s}.effort==efforts(i)))./sqrt(sum(bhvs{s}.effort==efforts(i)));
        bhvs{s}.dec_stake_mn(i) = mean(log(bhvs{s}.decision(bhvs{s}.stake==stakes(i))));
        bhvs{s}.dec_stake_se(i) = std(log(bhvs{s}.decision(bhvs{s}.stake==stakes(i))))./sqrt(sum(bhvs{s}.stake==stakes(i)));
    end
end

%% Full model including difficulty
% decision ease/difficulty did predict RTs, but check here reveals no
% effects on decision
fit_method = 'Laplace';
lme_full_str = 'decision_cur~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n)';
lme_full_noDEc_str = 'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n)';
lme_full_noDEp_str = 'decision_cur~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + (1|sbj_n)';
lme_full_noRc_str  = 'decision_cur~ effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n)';
lme_full_noESc_str = 'decision_cur~ reward_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n)';
lme_full = fitglme(good_tbl_prv.decision_cur,lme_full_str,'Distribution','binomial','FitMethod',fit_method);
lme_full_noDEc = fitglme(good_tbl_prv.decision_cur,lme_full_noDEc_str,'Distribution','binomial','FitMethod',fit_method);
lme_full_noDEp = fitglme(good_tbl_prv.decision_cur,lme_full_noDEp_str,'Distribution','binomial','FitMethod',fit_method);
lme_full_noRc  = fitglme(good_tbl_prv.decision_cur,lme_full_noRc_str,'Distribution','binomial','FitMethod',fit_method);
lme_full_noESc = fitglme(good_tbl_prv.decision_cur,lme_full_noESc_str,'Distribution','binomial','FitMethod',fit_method);

dec_dezc = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% NO improvement with decicsion ease
dec_dezp = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% NO, AIC worse with prv ease
dec_rewc = compare(lme_full_noRc,lme_full,'CheckNesting',true)
% YES, reward still very sig
dec_efSc = compare(lme_full_noESc,lme_full,'CheckNesting',true)
% YES, effortS still very sig

% Compute correlations between reward and decision ease
rew_rewchg_corr = nan(size(SBJs));
rew_rewchg_pval = nan(size(SBJs));
rew_ez_corr = nan(size(SBJs));
rew_ez_pval = nan(size(SBJs));
for s = 1:length(SBJs)
    % Comptue for reward and decision ease
    sbj_idx = good_tbl_all.decision_cur.sbj_n==s;
    [tmp,tmp_p] = corrcoef(good_tbl_all.decision_cur.reward_cur(sbj_idx),good_tbl_all.decision_cur.dec_ease_cur(sbj_idx));
    rew_ez_corr(s) = tmp(1,2);
    rew_ez_pval(s) = tmp_p(1,2);
    % Compute for reward (current) and reward change/contrast
    sbj_idx = good_tbl_grs.decision_cur.sbj_n==s;
    [tmp,tmp_p] = corrcoef(good_tbl_grs.decision_cur.reward_cur(sbj_idx),good_tbl_grs.decision_cur.reward_chg(sbj_idx));
    rew_rewchg_corr(s) = tmp(1,2);
    rew_rewchg_pval(s) = tmp_p(1,2);
end

%% Full model without difficulty
fit_method = 'Laplace';
model_formula_full = 'decision_cur ~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_noefSc  = 'decision_cur ~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_norewc  = 'decision_cur ~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_norewp  = 'decision_cur ~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_noefSp  = 'decision_cur ~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)';% + (1|trl_n_cur)';
dec_full = fitglme(good_tbl_prv.decision_cur,model_formula_full,'Distribution','binomial','FitMethod',fit_method);
dec_norc = fitglme(good_tbl_prv.decision_cur,model_formula_norewc,'Distribution','binomial','FitMethod',fit_method);
dec_noec = fitglme(good_tbl_prv.decision_cur,model_formula_noefSc,'Distribution','binomial','FitMethod',fit_method);
dec_norp = fitglme(good_tbl_prv.decision_cur,model_formula_norewp,'Distribution','binomial','FitMethod',fit_method);
dec_noep = fitglme(good_tbl_prv.decision_cur,model_formula_noefSp,'Distribution','binomial','FitMethod',fit_method);
dec_rew = compare(dec_norc,dec_full,'CheckNesting',true)
dec_efS = compare(dec_noec,dec_full,'CheckNesting',true)
dec_rewp = compare(dec_norp,dec_full,'CheckNesting',true)
dec_efSp = compare(dec_noep,dec_full,'CheckNesting',true)

model_formula_REcur = 'decision_cur ~ reward_cur + effortS_cur + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_noR  = 'decision_cur ~ effortS_cur + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_noE  = 'decision_cur ~ reward_cur + (1|sbj_n)';% + (1|trl_n_cur)';
dec_RE = fitglme(good_tbl_prv.decision_cur,model_formula_REcur,'Distribution','binomial','FitMethod',fit_method);
dec_noR  = fitglme(good_tbl_prv.decision_cur,model_formula_noR,'Distribution','binomial','FitMethod',fit_method);
dec_noE  = fitglme(good_tbl_prv.decision_cur,model_formula_noE,'Distribution','binomial','FitMethod',fit_method);
dec_R = compare(dec_noR,dec_RE,'CheckNesting',true)
dec_E = compare(dec_noE,dec_RE,'CheckNesting',true)

%% Test previous reward interactions with reward/effort
fit_method = 'Laplace';
model_formula_full = 'decision_cur ~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
dec_full = fitglme(good_tbl_prv.decision_cur,model_formula_full,'Distribution','binomial','FitMethod',fit_method);

lme_full_pRcRint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
lme_full_pRcEint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
lme_full_pRpEint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
dec_pRcRint_p = compare(dec_full,lme_full_pRcRint,'CheckNesting',true)
dec_pRcEint_p = compare(dec_full,lme_full_pRcEint,'CheckNesting',true)
dec_pRpEint_p = compare(dec_full,lme_full_pRpEint,'CheckNesting',true)
% dec_pRcR_vs_pRcE_p = compare(lme_full_pRcEint,lme_full_pRcRint,'NSim',1000)

lme_full_pRcREint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur*reward_prv + effortS_cur*reward_prv + effortS_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
dec_pRcREint_p = compare(lme_full_pRcRint,lme_full_pRcREint,'CheckNesting',true)
dec_pRcREint_p2 = compare(lme_full_pRcEint,lme_full_pRcREint,'CheckNesting',true)

% Post-hoc tests for differences on median splits
rew_prv_lo_idx = good_tbl_prv.decision_cur.reward_prv<-0.5;
rew_prv_hi_idx = good_tbl_prv.decision_cur.reward_prv>0.5;
model_formula_nocR = 'decision_cur ~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';
model_formula_nocE = 'decision_cur ~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)';
dec_pR_lo = fitglme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),model_formula_full,'Distribution','binomial','FitMethod',fit_method);
dec_pR_hi = fitglme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),model_formula_full,'Distribution','binomial','FitMethod',fit_method);
dec_pR_lo_nocR = fitglme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),model_formula_nocR,'Distribution','binomial','FitMethod',fit_method);
dec_pR_hi_nocR = fitglme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),model_formula_nocR,'Distribution','binomial','FitMethod',fit_method);
dec_pR_lo_nocE = fitglme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),model_formula_nocE,'Distribution','binomial','FitMethod',fit_method);
dec_pR_hi_nocE = fitglme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),model_formula_nocE,'Distribution','binomial','FitMethod',fit_method);
dec_pRlo_cR = compare(dec_pR_lo_nocR,dec_pR_lo,'CheckNesting',true)
dec_pRhi_cR = compare(dec_pR_hi_nocR,dec_pR_hi,'CheckNesting',true)
dec_pRlo_cE = compare(dec_pR_lo_nocE,dec_pR_lo,'CheckNesting',true)
dec_pRhi_cE = compare(dec_pR_hi_nocE,dec_pR_hi,'CheckNesting',true)


% Test if PFC theta is better than previous reward
% lme_full_TcRint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:PFC_theta + (1|sbj_n)',...
%     'Distribution','binomial','FitMethod',fit_method);
% lme_full_TcEint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:PFC_theta + (1|sbj_n)',...
%     'Distribution','binomial','FitMethod',fit_method);
% lme_full_TpEint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:PFC_theta + (1|sbj_n)',...
%     'Distribution','binomial','FitMethod',fit_method);
% dec_TcRint_p = compare(dec_full,lme_full_TcRint,'CheckNesting',true)
% dec_TcEint_p = compare(dec_full,lme_full_TcEint,'CheckNesting',true)
% dec_TpEint_p = compare(dec_full,lme_full_TpEint,'CheckNesting',true)
% All no, but the previous effort one is actually close (p=0.053)...

% Plot partial dependence
fn_plot_LMM_interaction_partial_dependence(lme_full_pRcRint,'reward_cur','reward_prv');
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_interaction_partial_dependence(lme_full_pRcEint,'effortS_cur','reward_prv');
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
% fn_plot_LMM_interaction_partial_dependence(lme_full_pRpEint,'effortS_prv','reward_prv');
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot partial dependence in model with both interactions
% fn_plot_LMM_interaction_partial_dependence(lme_full_pRcREint,'reward_cur','reward_prv');
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '_pRcREint.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
% fn_plot_LMM_interaction_partial_dependence(lme_full_pRcREint,'effortS_cur','reward_prv');
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '_pRcREint.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot tertile line plots
% fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','reward_cur','decision_cur',3,5);
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
% fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_prv','decision_cur',3,5);
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
% fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','decision_cur',3,5);
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot median splits as line plots 
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','reward_cur','decision_cur',2,5);
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','decision_cur',2,5);
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_prv','decision_cur',2,5);

% Plot median split bar plots
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.decision_cur,'reward_prv','reward_cur','decision_cur');
ylabel('Offers Accepted (%)');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_gratton_diff_lines_sbj(good_tbl_prv.decision_cur,'reward_prv','reward_cur','decision_cur');
ylabel('Offers Accepted (%) Hi-Lo Cur Reward');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','decision_cur');
ylabel('Offers Accepted (%)');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_gratton_diff_lines_sbj(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','decision_cur');
ylabel('Offers Accepted (%) Hi-Lo Cur EffortS');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Test reward:effort interactions
formula_cREint = 'decision_cur ~ reward_cur + effortS_cur + reward_cur:effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';
dec_full = fitglme(good_tbl_prv.decision_cur,model_formula_full,'Distribution','binomial','FitMethod',fit_method);
lme_cREint = fitglme(good_tbl_prv.decision_cur,formula_cREint,'Distribution','binomial','FitMethod',fit_method);
dec_cREint_p = compare(dec_full,lme_cREint,'CheckNesting',true)
% Nope, don't add R:E interaction

%% pAccept previous reward interactions
fit_method = 'ML';
pacc_formula_RE_full = 'pAccept_cur ~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
pacc_full = fitlme(good_tbl_prv.decision_cur,pacc_formula_RE_full,'FitMethod',fit_method);
pacc_lme_nocR = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)','FitMethod',fit_method);
pacc_lme_nocE = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)','FitMethod',fit_method);
pacc_lme_nopR = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)','FitMethod',fit_method);
pacc_lme_nopE = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)','FitMethod',fit_method);
pacc_lme_full_pRcRint = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)',...
    'FitMethod',fit_method);
pacc_lme_full_pRcEint = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:reward_prv + (1|sbj_n)',...
    'FitMethod',fit_method);
pacc_lme_full_pRpEint = fitlme(good_tbl_prv.decision_cur,'pAccept_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:reward_prv + (1|sbj_n)',...
    'FitMethod',fit_method);
pacc_cR_p = compare(pacc_lme_nocR,pacc_full,'CheckNesting',true)
pacc_cE_p = compare(pacc_lme_nocE,pacc_full,'CheckNesting',true)
pacc_pR_p = compare(pacc_lme_nopR,pacc_full,'CheckNesting',true)
pacc_pE_p = compare(pacc_lme_nopE,pacc_full,'CheckNesting',true)
pacc_pRcRint_p = compare(pacc_full,pacc_lme_full_pRcRint,'CheckNesting',true)
pacc_pRcEint_p = compare(pacc_full,pacc_lme_full_pRcEint,'CheckNesting',true)
pacc_pRpEint_p = compare(pacc_full,pacc_lme_full_pRpEint,'CheckNesting',true)

% Post-hoc tests for differences on median splits
pacc_formula_nocR = 'pAccept_cur ~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';
pacc_formula_nocE = 'pAccept_cur ~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)';
pacc_pR_lo = fitlme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),pacc_formula_RE_full,'FitMethod',fit_method);
pacc_pR_hi = fitlme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),pacc_formula_RE_full,'FitMethod',fit_method);
pacc_pR_lo_nocR = fitlme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),pacc_formula_nocR,'FitMethod',fit_method);
pacc_pR_hi_nocR = fitlme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),pacc_formula_nocR,'FitMethod',fit_method);
pacc_pR_lo_nocE = fitlme(good_tbl_prv.decision_cur(rew_prv_lo_idx,:),pacc_formula_nocE,'FitMethod',fit_method);
pacc_pR_hi_nocE = fitlme(good_tbl_prv.decision_cur(rew_prv_hi_idx,:),pacc_formula_nocE,'FitMethod',fit_method);
pacc_pRlo_cR = compare(pacc_pR_lo_nocR,pacc_pR_lo,'CheckNesting',true)
pacc_pRhi_cR = compare(pacc_pR_hi_nocR,pacc_pR_hi,'CheckNesting',true)
pacc_pRlo_cE = compare(pacc_pR_lo_nocE,pacc_pR_lo,'CheckNesting',true)
pacc_pRhi_cE = compare(pacc_pR_hi_nocE,pacc_pR_hi,'CheckNesting',true)

% Plot partial dependence
fn_plot_LMM_interaction_partial_dependence(pacc_lme_full_pRcRint,'reward_cur','reward_prv');
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_interaction_partial_dependence(pacc_lme_full_pRcEint,'effortS_cur','reward_prv');
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot median splits as line plots 
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','reward_cur','pAccept_cur',2,5);
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','pAccept_cur',2,5);
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.decision_cur,'reward_prv','effortS_prv','pAccept_cur',2,5);
set(gca,'FontSize',20);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot median split bar plots
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.decision_cur,'reward_prv','reward_cur','pAccept_cur');
ylabel('Model p(Accept)'); ylim([0 1]);
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_gratton_diff_lines_sbj(good_tbl_prv.decision_cur,'reward_prv','reward_cur','pAccept_cur');
ylabel('Model p(Accept) Hi-Lo Cur Reward');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','pAccept_cur');
ylabel('Model p(Accept)'); ylim([0 1]);
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_gratton_diff_lines_sbj(good_tbl_prv.decision_cur,'reward_prv','effortS_cur','pAccept_cur');
ylabel('Model p(Accept) Hi-Lo Cur EffS');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

%% Compare SV vs. reward + effort
fit_method = 'Laplace';
model_formula_reweffS = 'decision_cur ~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_sv  = 'decision_cur ~ SV_cur + SV_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_nosv   = 'decision_cur ~ SV_prv + (1|sbj_n)';% + (1|trl_n_cur)';
model_formula_nosvp  = 'decision_cur ~ SV_cur + (1|sbj_n)';% + (1|trl_n_cur)';
dec_reweffS = fitglme(good_tbl_prv.decision_cur,model_formula_reweffS,'Distribution','binomial','FitMethod',fit_method);
dec_sv      = fitglme(good_tbl_prv.decision_cur,model_formula_sv,'Distribution','binomial','FitMethod',fit_method);
dec_nosv    = fitglme(good_tbl_prv.decision_cur,model_formula_nosv,'Distribution','binomial','FitMethod',fit_method);
dec_nosvp   = fitglme(good_tbl_prv.decision_cur,model_formula_nosvp,'Distribution','binomial','FitMethod',fit_method);
dec_reweff_vs_sv = compare(dec_sv,dec_reweffS)
dec_SV_pval = compare(dec_nosv,dec_sv,'CheckNesting',true)
dec_SVp_pval = compare(dec_nosvp,dec_sv,'CheckNesting',true)

%% Test previous reward interactions with subjective value
fit_method = 'Laplace';
model_formula_sv  = 'decision_cur ~ SV_cur + SV_prv + (1|sbj_n)';% + (1|trl_n_cur)';
dec_sv_full = fitglme(good_tbl_prv.decision_cur,model_formula_sv,'Distribution','binomial','FitMethod',fit_method);
lme_full_SVpR = fitglme(good_tbl_prv.decision_cur,'decision_cur~ SV_cur + SV_prv + reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
lme_full_pRSVint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ SV_cur + SV_prv + SV_cur:reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
lme_full_pRpSVint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ SV_cur + SV_prv + SV_prv:reward_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
lme_full_pSVSVint = fitglme(good_tbl_prv.decision_cur,'decision_cur~ SV_cur*SV_prv + (1|sbj_n)',...
    'Distribution','binomial','FitMethod',fit_method);
dec_SVpR_p = compare(dec_sv_full,lme_full_SVpR,'CheckNesting',true)
% Nope
dec_pRSVint_p = compare(dec_sv_full,lme_full_pRSVint,'CheckNesting',true)
% Nope
dec_pRpSVint_p = compare(dec_sv_full,lme_full_pRpSVint,'CheckNesting',true)
% Nope
dec_pSVSVint_p = compare(dec_sv_full,lme_full_pSVSVint,'CheckNesting',true)
% Nope

%% Test decision ~ current individual task features
fit_method = 'Laplace';
% Decision ~ reward_cur
lme0 = fitglme(good_tbl_all.decision_cur,'decision_cur~ 1 + (1|sbj_n)','Distribution','binomial','FitMethod',fit_method);
lme_rewc = fitglme(good_tbl_all.decision_cur,'decision_cur~ reward_cur + (1|sbj_n)','Distribution','binomial','FitMethod',fit_method);
% dec_rew = compare(lme0,lme_rewc,'CheckNesting',true)%,'NSim',1000)

fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.decision_cur,'reward_cur','decision_cur',...
    lme_rewc,lme_rewc.Coefficients.pValue(2),n_quantiles,'logistic_fit',1);
xlabel('Reward (z)');
ylabel('Decision');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

lme_dezc = fitglme(good_tbl_all.decision_cur,'decision_cur~ dec_ease_cur + (1|sbj_n)','Distribution','binomial','FitMethod',fit_method);
fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.decision_cur,'dec_ease_cur','decision_cur',...
    lme_dezc,lme_dezc.Coefficients.pValue(2),n_quantiles);%,'logistic_fit',1);
xlabel('Reward (z)');
ylabel('Decision');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Test GRS (not at all close)
% lme0    = fitglme(good_tbl_grs.decision_cur,'decision_cur~ 1 + (1|sbj_n)','Distribution','binomial','FitMethod',fit_method);
% lme_grs = fitglme(good_tbl_grs.decision_cur,'decision_cur~ grs + (1|sbj_n)','Distribution','binomial','FitMethod',fit_method);
% dec_grs = compare(lme0,lme_grs,'CheckNesting',true)%,'NSim',1000)

%% Neural power predict decision
% lme_fullpow = fitglme(good_tbl_all.decision_cur,...
%     'decision_cur~ PFC_theta + BG_theta + PFC_betalo + BG_betalo + (1|sbj_n)','Distribution','binomial');
% BG_betalo marginal (p=0.06), nothing else close
