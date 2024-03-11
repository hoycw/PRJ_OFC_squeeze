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
scat_sz = 40;

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);

fig_dir   = [prj_dir 'results/bhv/LMM/' an_id '/' stat_id '/lRT/'];
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
rt_vars = {'rt_cur','logrt_cur','rt_prv','logrt_prv'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(rt_vars)
    % Identify outliers
    out_idx_all.(rt_vars{f}) = abs(table_all.(rt_vars{f}))>st.outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(rt_vars{f}) = table_all(~out_idx_all.(rt_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',rt_vars{f});
    if any(out_idx_all.(rt_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(rt_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(rt_vars{f})),rt_vars{f});
        fprintf(2,'%.2f, ',table_all.(rt_vars{f})(out_idx_all.(rt_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',rt_vars{f},st.outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        rt_vars{f},size(good_tbl_all.(rt_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create previous trial and GRS tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
tbl_fields = good_tbl_all.(rt_vars{1}).Properties.VariableNames;
for p = 1:length(rt_vars)
    prv_nan_idx = isnan(good_tbl_prv.(rt_vars{p}).SV_prv);
    good_tbl_prv.(rt_vars{p})(prv_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_prv.(rt_vars{p}).(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
            error(['NaN is table_prv.' tbl_fields{f}]);
        end
    end
end

% Toss NaNs from grs table
good_tbl_grs = good_tbl_all;
for p = 1:length(rt_vars)
    grs_nan_idx = isnan(good_tbl_grs.(rt_vars{p}).grs);
    good_tbl_grs.(rt_vars{p})(grs_nan_idx,:) = [];
    for f = 1:length(tbl_fields)
        if any(isnan(good_tbl_grs.(rt_vars{p}).(tbl_fields{f})))
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
    bhvs{s}.lrt_stake_mn  = nan([5 1]);
    bhvs{s}.lrt_stake_se  = nan([5 1]);
%     bhvs{s}.lrt_effort_mn = nan([5 1]);
%     bhvs{s}.lrt_effort_se = nan([5 1]);
    for i = 1:5
        bhvs{s}.SV_effort_mn(i) = mean(bhvs{s}.SV(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_effort_se(i) = std(bhvs{s}.SV(bhvs{s}.effort==efforts(i)))./sqrt(sum(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_stake_mn(i) = mean(bhvs{s}.SV(bhvs{s}.stake==stakes(i)));
        bhvs{s}.SV_stake_se(i) = std(bhvs{s}.SV(bhvs{s}.stake==stakes(i)))./sqrt(sum(bhvs{s}.stake==stakes(i)));
%         bhvs{s}.lrt_effort_mn(i) = mean(bhvs{s}.rt(bhvs{s}.effort==efforts(i)));
%         bhvs{s}.lrt_effort_se(i) = std(bhvs{s}.rt(bhvs{s}.effort==efforts(i)))./sqrt(sum(bhvs{s}.effort==efforts(i)));
        bhvs{s}.lrt_stake_mn(i) = mean(log(bhvs{s}.rt(bhvs{s}.stake==stakes(i))));
        bhvs{s}.lrt_stake_se(i) = std(log(bhvs{s}.rt(bhvs{s}.stake==stakes(i))))./sqrt(sum(bhvs{s}.stake==stakes(i)));
    end
end

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewc = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewp = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)');
lme_full_noeffc = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_noeffp = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');

logrt_cur_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
logrt_cur_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
logrt_cur_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
logrt_cur_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

% Sup. Fig. 1c: Plot logRT ~ reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.logrt_cur,'reward_cur','logrt_cur',...
    lme_full,logrt_cur_rewc.pValue(2),n_quantiles);
xlabel('Reward (z)');
ylabel('log(RT)');
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

%% Test previous reward interactions
lme_full = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_pRcRint = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)');
lme_full_pRcEint = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:reward_prv + (1|sbj_n)');
lme_full_pRpEint = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:reward_prv + (1|sbj_n)');
logrt_cur_pRcRint_p = compare(lme_full,lme_full_pRcRint,'CheckNesting',true) % no
logrt_cur_pRcEint_p = compare(lme_full,lme_full_pRcEint,'CheckNesting',true) % no
logrt_cur_pRpEint_p = compare(lme_full,lme_full_pRpEint,'CheckNesting',true) % no

%% Test effect of decision ease
lme_ez        = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ dec_ease_cur + dec_ease_prv + (1|sbj_n)');
lme_ez_noDEc  = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ dec_ease_prv + (1|sbj_n)');
lme_ez_noDEp  = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ dec_ease_cur + (1|sbj_n)');

lrt_ezc = compare(lme_ez_noDEc,lme_ez,'CheckNesting',true) % yes, significant
lrt_ezp = compare(lme_ez_noDEp,lme_ez,'CheckNesting',true) % yes, significant

%% Compare effect of reward and decision ease in same model
lme_fullez        = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effort_cur + dec_ease_cur + reward_prv + effort_prv + dec_ease_prv + (1|sbj_n)');
lme_fullez_noDEc  = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effort_cur + reward_prv + effort_prv + dec_ease_prv + (1|sbj_n)');
lme_fullez_noDEp  = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effort_cur + dec_ease_cur + reward_prv + effort_prv + (1|sbj_n)');
lme_fullez_norewc = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ effort_cur + dec_ease_cur + reward_prv + effort_prv + dec_ease_prv + (1|sbj_n)');

lrt_dec = compare(lme_fullez_noDEc,lme_fullez,'CheckNesting',true) % yes, ease still significant
lrt_dep = compare(lme_fullez_noDEp,lme_fullez,'CheckNesting',true) % yes, previous ease still significant
lrt_rew = compare(lme_fullez_norewc,lme_fullez,'CheckNesting',true) % no, reward is no longer significant

% Sup. Fig. 1b: gratton bar plot for decision ease
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.logrt_cur,'dec_ease_prv','dec_ease_cur','logrt_cur');
ylabel('log(RT)');
set(gca,'FontSize',18);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Compute correlations between reward and decision ease
rew_ez_corr = nan(size(SBJs));
rew_ez_pval = nan(size(SBJs));
for s = 1:length(SBJs)
    sbj_idx = good_tbl_all.logrt_cur.sbj_n==s;
    [tmp,tmp_p] = corrcoef(good_tbl_all.logrt_cur.reward_cur(sbj_idx),good_tbl_all.logrt_cur.dec_ease_cur(sbj_idx));
    rew_ez_corr(s) = tmp(1,2);
    rew_ez_pval(s) = tmp_p(1,2);
end

% Confirm ease overrules reward with single predictors
% lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
% lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ dec_ease_cur + (1|sbj_n)');
% lme3 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + dec_ease_cur + (1|sbj_n)');
% lme4 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + effort_cur + dec_ease_cur + (1|sbj_n)');
% lrt_rew = compare(lme0,lme1) % yes
% lrt_ease = compare(lme0,lme2) % yes
% lrt_rewez = compare(lme1,lme3) % yes
% lrt_ezrew = compare(lme2,lme3) % no

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ SV_cur + SV_prv + (1|sbj_n)');
logrt_cur_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ SV_cur + (1|sbj_n)');
lme3 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + SV_cur + (1|sbj_n)');    % neither significant when both in model
lrt_rew = compare(lme0,lme1)%,'NSim',1000)
lrt_svc = compare(lme0,lme2)
lrt_rewsv = compare(lme1,lme3,'NSim',1000) % adding SV to reward does not improve model fit
lrt_rew_vs_sv = compare(lme1,lme2,'NSim',1000) % no significant difference in model fit between rew_cur and SV_cur

%% Test log(RT) vs. RT modeling
% Reward model (this result holds for effort and SV too)
lme0 = fitlme(good_tbl_all.rt_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.rt_cur,'rt_cur~ reward_cur + (1|sbj_n)');
rt_rew = compare(lme0,lme1)%,'NSim',1000)
lme0_log = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1_log = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
logrt_rew = compare(lme0_log,lme1_log)%,'NSim',1000)

fig_name = 'GRP_RT_log_transform_QA';
figure('Name',fig_name);
subplot(2,2,1);
histogram(bhvs{s}.rt,'FaceColor','b');
title('RTs'); set(gca,'FontSize',16);
subplot(2,2,2);
histogram(log(bhvs{s}.rt),'FaceColor','r');
title('log(RTs)');  set(gca,'FontSize',16);
subplot(2,2,3);
qqplot(lme1.residuals); set(gca,'FontSize',16);
title('QQ Plot: RT vs. Reward');
subplot(2,2,4);
qqplot(lme1_log.residuals); set(gca,'FontSize',16);
title('QQ Plot: log(RT) vs. Reward');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Salience models
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ dec_ease_cur + (1|sbj_n)');
rt_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_dec_ease = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% ease is better than abs(SV) according to AIC and BIC

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Decision model
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ decision_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ pAccept_cur + (1|sbj_n)');
rt_dec = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_pAcc = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

%% Beta
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ PFC_betalo + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ BG_betalo + (1|sbj_n)');
lrt_pfcbeta = compare(lme0,lme1)%,'NSim',1000)
lrt_bgbeta = compare(lme0,lme2)%,'NSim',1000)

