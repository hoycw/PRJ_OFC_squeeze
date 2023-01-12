%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_madS8t0_f2t40';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
outlier_thresh = 4;
n_quantiles = 5;

if contains(an_id,'_S')
    an_lim = [0.5 1.5];
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end

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

win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');
table_name = [an_id win_str norm_bhv_str norm_nrl_str];
fig_dir   = [prj_dir 'results/bhv/LMM/' table_name out_thresh_str '/lRT/'];
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
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Toss outliers
rt_vars = {'rt_cur','logrt_cur','rt_prv','logrt_prv'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(rt_vars)
    % Identify outliers
    out_idx_all.(rt_vars{f}) = abs(table_all.(rt_vars{f}))>outlier_thresh;
    
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
        fprintf('No bad trials for %s with threshold %d\n',rt_vars{f},outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        rt_vars{f},size(good_tbl_all.(rt_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create current and previous trial tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
for p = 1:length(rt_vars)
    prv_nan_idx = isnan(good_tbl_prv.(rt_vars{p}).SV_prv);
    good_tbl_prv.(rt_vars{p})(prv_nan_idx,:) = [];
    prv_fields = good_tbl_prv.(rt_vars{p}).Properties.VariableNames;
    for f = 1:length(prv_fields)
        if any(isnan(good_tbl_prv.(rt_vars{p}).(prv_fields{f}))); error(['NaN is good_tbl_prv.' prv_fields{f}]); end
    end
end

%% Behavioral model fitting
%   SV_physical(t) = R(t)?k*E(t)^2
%   softmax: exp(B * Q1) / exp(B * Q1) + exp(B * Q2)
%   here, Q1 is SV and Q2 is nothing (reject choice, no alternative, so just B term)
% decisionfun=@(p) norm( (exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))) ./ ...
%     (exp(p(1)) + exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))))) - bhv.decision);
% [par, fit]=fminsearch(decisionfun, [1,1]);
% 
% SV_fn    = @(k) bhv.stake-(k*(bhv.effort).^2);
% EFF_fn   = @(k) (k*(bhv.effort).^2);
% bhv.SV   = SV_fn(par(2));
% bhv.EFFs = EFF_fn(par(2));
% bhv.p_accept = (exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2))) ./...
%     (exp(par(1)) + exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2)))));

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

%% Test decision ~ task features
% dec0 = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ 1 + (1|sbj_n)','Distribution','binomial');
% decr = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ reward_cur + (1|sbj_n)','Distribution','binomial');
% dec_rew = compare(dec0,decr,'CheckNesting',true)%,'NSim',1000)
% dece  = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ effort_cur + (1|sbj_n)','Distribution','binomial');
% deceS = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ effortS_cur + (1|sbj_n)','Distribution','binomial');
% dec_eff = compare(dec0,dece,'CheckNesting',true)%,'NSim',1000)
% dec_effS = compare(dec0,deceS,'CheckNesting',true)%,'NSim',1000)
% % dec_effS = compare(deceS,decreS,'CheckNesting',true)%,'NSim',1000)
% %   regular effort seems to be better predictor than subjective effort
% decsv  = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ SV_cur + (1|sbj_n)','Distribution','binomial');
% 
% decre  = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ effort_cur + reward_cur + (1|sbj_n)','Distribution','binomial');
% decreS = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ effortS_cur + reward_cur + (1|sbj_n)','Distribution','binomial');
% dec_re_vs_sv = compare(decsv,decre)
% 
% decrei = fitglme(good_tbl_all.logrt_cur,'decision_cur ~ effort_cur*reward_cur + (1|sbj_n)','Distribution','binomial');
% dec_effS = compare(dec0,deceS,'CheckNesting',true)%,'NSim',1000)
% 
% % Full model
% model_formula_full = 'decision_cur ~ effort_cur + reward_cur + dec_diff_cur + (1|sbj_n)';
% model_formula_noE  = 'decision_cur ~ reward_cur + dec_diff_cur + (1|sbj_n)';
% model_formula_noR  = 'decision_cur ~ effort_cur + dec_diff_cur + (1|sbj_n)';
% model_formula_noD  = 'decision_cur ~ reward_cur + effort_cur + (1|sbj_n)';
% decred = fitglme(good_tbl_all.logrt_cur,model_formula_full,'Distribution','binomial');
% decrd = fitglme(good_tbl_all.logrt_cur,model_formula_noE,'Distribution','binomial');
% deced = fitglme(good_tbl_all.logrt_cur,model_formula_noR,'Distribution','binomial');
% decre = fitglme(good_tbl_all.logrt_cur,model_formula_noD,'Distribution','binomial');
% dec_rew = compare(deced,decred,'CheckNesting',true)
% dec_eff = compare(decrd,decred,'CheckNesting',true)

%% Full model of everything
lme0 = fitlme(good_tbl_all.rt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.rt_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.rt_cur,'logrt_cur~ dec_diff_cur + (1|sbj_n)');
lme3 = fitlme(good_tbl_all.rt_cur,'logrt_cur~ reward_cur + dec_diff_cur + (1|sbj_n)');
lme4 = fitlme(good_tbl_all.rt_cur,'logrt_cur~ reward_cur + effort_cur + dec_diff_cur + (1|sbj_n)');
lrt_rew = compare(lme0,lme1)%,'NSim',1000)
lrt_ease = compare(lme0,lme2)
lrt_rewez = compare(lme1,lme3) % adding decision ease to reward improves model fit
lrt_ezrew = compare(lme2,lme3) % adding reward to decision ease is not significantly better

lme_full = fitlme(good_tbl_prv.rt_cur,'logrt_cur~ reward_cur + effort_cur + dec_diff_cur + reward_prv + effort_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noDEc = fitlme(good_tbl_prv.rt_cur,'logrt_cur~ reward_cur + effort_cur + reward_prv + effort_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noDEp = fitlme(good_tbl_prv.rt_cur,'logrt_cur~ reward_cur + effort_cur + dec_diff_cur + reward_prv + effort_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewc = fitlme(good_tbl_prv.rt_cur,'logrt_cur~ effort_cur + dec_diff_cur + reward_prv + effort_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');

lrt_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true) % adding reward to decision ease is not significantly better
lrt_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true) % adding reward to decision ease is not significantly better
lrt_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true) % adding reward to decision ease is not significantly better

% Compute correlations between reward and decision ease
rew_ez_corr = nan(size(SBJs));
rew_ez_pval = nan(size(SBJs));
for s = 1:length(SBJs)
    sbj_idx = good_tbl_all.rt_cur.sbj_n==s;
    [tmp,tmp_p] = corrcoef(good_tbl_all.rt_cur.reward_cur(sbj_idx),good_tbl_all.rt_cur.dec_diff_cur(sbj_idx));
    rew_ez_corr(s) = tmp(1,2);
    rew_ez_pval(s) = tmp_p(1,2);
end

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

%% Current Reward, Effort,a nd SV models
scat_sz = 40;

% Reward model
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
%     AIC       BIC       LogLikelihood    Deviance
%     1631.2    1648.6    -811.61          1623.2  
rt_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% Plot logRT ~ reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_all.logrt_cur,'reward_cur','logrt_cur',...
    lme1,rt_rew.pValue(2),n_quantiles);
xlabel('Reward (z)');
ylabel('log RT (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

fig_name = 'GRP_lRT_LMM_results';
figure('Name',fig_name);
subplot(2,2,1); hold on;
for s = 1:length(SBJs)
%     scatter(good_tbl_all.logrt_cur.reward_cur(good_tbl_all.logrt_cur.sbj_n==s),good_tbl_all.logrt_cur.logrt_cur(good_tbl_all.logrt_cur.sbj_n==s),...
%         scat_sz,sbj_colors(s,:));
    errorbar(zscore(stakes),zscore(bhvs{s}.lrt_stake_mn),bhvs{s}.lrt_stake_se,'Color',sbj_colors(s,:));
end
xvals = min(good_tbl_all.logrt_cur.reward_cur):0.01:max(good_tbl_all.logrt_cur.reward_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Reward (z)');
xlim([-2 2]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_rew.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

% logRTs ~ Effort current
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ effort_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ effortS_cur + (1|sbj_n)');
% rt_effort_EFF = compare(lme1,lme2)%,'NSim',1000)
%   effort is better predictor than EFFs (subjective)
rt_effort = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Subjective Value model
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ SV_cur + (1|sbj_n)');
%     AIC       BIC     LogLikelihood    Deviance
%     1631.6    1649    -811.8           1623.6  
% lme3 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + effort_cur + (1|sbj_n)');
%     AIC       BIC       LogLikelihood    Deviance
%     1631.5    1653.3    -810.76          1621.5  
% lme4 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur*effortS_cur + (1|sbj_n)');
%     AIC       BIC       LogLikelihood    Deviance
%     1629.7    1655.8    -808.86          1617.7  
% lme42 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur:effortS_cur + (1|sbj_n)');
%     AIC       BIC       LogLikelihood    Deviance
%     1631.1    1648.5    -811.53          1623.1  
% lme5 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ absSV_cur + (1|sbj_n)');
%     AIC       BIC     LogLikelihood    Deviance
%     1609.6    1627    -800.79          1601.6  
% lme6 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ reward_cur + absSV_cur + (1|sbj_n)');
%     AIC       BIC       LogLikelihood    Deviance
%     1611.5    1633.3    -800.75          1601.5  
rt_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
subplot(2,2,2); hold on;
for s = 1:length(SBJs)
    scatter(good_tbl_all.logrt_cur.SV_cur(good_tbl_all.logrt_cur.sbj_n==s),good_tbl_all.logrt_cur.logrt_cur(good_tbl_all.logrt_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(good_tbl_all.logrt_cur.SV_cur):0.01:max(good_tbl_all.logrt_cur.SV_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Subjective Value (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_sv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

%% Salience models
lme0 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.logrt_cur,'logrt_cur~ dec_diff_cur + (1|sbj_n)');
rt_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_dec_diff = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
lrt_best_diff = compare(lme1,lme2,'NSim',1000)
subplot(2,2,3); hold on;
for s = 1:length(SBJs)
    scatter(good_tbl_all.logrt_cur.absSV_cur(good_tbl_all.logrt_cur.sbj_n==s),good_tbl_all.logrt_cur.logrt_cur(good_tbl_all.logrt_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(good_tbl_all.logrt_cur.absSV_cur):0.01:max(good_tbl_all.logrt_cur.absSV_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('absolute SV (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_abssv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);
subplot(2,2,4); hold on;
for s = 1:length(SBJs)
    scatter(good_tbl_all.logrt_cur.dec_diff_cur(good_tbl_all.logrt_cur.sbj_n==s),good_tbl_all.logrt_cur.logrt_cur(good_tbl_all.logrt_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(good_tbl_all.logrt_cur.dec_diff_cur):0.01:max(good_tbl_all.logrt_cur.dec_diff_cur);
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Decision Ease (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(rt_dec_diff.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

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

%% Previous trial predictors
% Previous Reward model
lme0 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ reward_prv + (1|sbj_n)');
rt_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Previous Effort model
lme0 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ effort_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ effortS_prv + (1|sbj_n)');
rt_effort_EFF_prv = compare(lme1,lme2)%,'NSim',1000)
rt_effort_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_effortS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Previous Subjective Value model
lme0 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ SV_prv + (1|sbj_n)');
rt_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Previous salience model
lme0 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ dec_diff_prv + (1|sbj_n)');
lme3 = fitlme(good_tbl_prv.logrt_cur,'logrt_cur~ logrt_prv + (1|sbj_n)');
rt_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_dec_diff_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

fig_name = 'GRP_lRT_LMM_results_prv_salience';
figure('Name',fig_name);
subplot(1,2,1); hold on;
for s = 1:length(SBJs)
    scatter(good_tbl_prv.logrt_cur.absSV_prv(good_tbl_prv.logrt_cur.sbj_n==s),good_tbl_prv.logrt_cur.logrt_cur(good_tbl_prv.logrt_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(good_tbl_prv.logrt_cur.absSV_prv):0.01:max(good_tbl_prv.logrt_cur.absSV_prv);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Previous abs(SV) (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_abssv_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);
subplot(1,2,2); hold on;
for s = 1:length(SBJs)
    scatter(good_tbl_prv.logrt_cur.dec_diff_prv(good_tbl_prv.logrt_cur.sbj_n==s),good_tbl_prv.logrt_cur.logrt_cur(good_tbl_prv.logrt_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
    
end
xvals = min(good_tbl_prv.logrt_cur.dec_diff_prv):0.01:max(good_tbl_prv.logrt_cur.dec_diff_prv);
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Previous Decision Ease (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(rt_dec_diff_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Gratton-style line plot
fig_name = 'GRP_lRT_LMM_results_logRT_salience_gratton';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.6 0.6]);
subplot(1,2,1); hold on;
prv_low_idx = good_tbl_prv.logrt_cur.absSV_prv<median(good_tbl_prv.logrt_cur.absSV_prv);
cur_low_idx = good_tbl_prv.logrt_cur.absSV_cur<median(good_tbl_prv.logrt_cur.absSV_cur);

lL_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & ~cur_low_idx));
lL_sem = std(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('log RT (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. abs(SV)','High Prev. abs(SV)'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. abs(SV)','High Curr. abs(SV)'},'Location','northwest');
title('Effects of Previous/Current abs(SV)');
set(gca,'FontSize',16);

subplot(1,2,2); hold on;
prv_low_idx = good_tbl_prv.logrt_cur.dec_diff_prv<median(good_tbl_prv.logrt_cur.dec_diff_prv);
cur_low_idx = good_tbl_prv.logrt_cur.dec_diff_cur<median(good_tbl_prv.logrt_cur.dec_diff_cur);

lL_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & ~cur_low_idx));
lL_sem = std(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.logrt_cur.logrt_cur(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.logrt_cur.logrt_cur(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('log RT (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. Ease','High Prev. Ease'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. Ease','High Curr. Ease'},'Location','northwest');
title('Effects of Previous/Current Decision Ease');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Does previous trial subjective value predict anything useful?
% error('!!! need to check/toss these variable outliers');
% % Current subjective value
% lme0 = fitlme(table_all,'SV_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'SV_cur~ SV_prv + (1|sbj_n)');
% sv_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Current probability of accept
% lme0 = fitlme(table_all,'pAccept_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'pAccept_cur~ SV_prv + (1|sbj_n)');
% pAcc_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Current decision
% lme0 = fitlme(table_all,'decision_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'decision_cur~ SV_prv + (1|sbj_n)');
% dec_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Current salience
% lme0 = fitlme(table_all,'absSV_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'absSV_cur~ SV_prv + (1|sbj_n)');
% abssv_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(table_all,'dec_diff_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'dec_diff_cur~ SV_prv + (1|sbj_n)');
% dec_diff_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Current reward
% lme0 = fitlme(table_all,'reward_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'reward_cur~ SV_prv + (1|sbj_n)');
% rew_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Current Effort
% lme0 = fitlme(table_all,'effort_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'effort_cur~ SV_prv + (1|sbj_n)');
% effort_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(table_all,'effortS_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'effortS_cur~ SV_prv + (1|sbj_n)');
% effortS_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % Why? check autocorrelation of effortS
% lme0 = fitlme(table_all,'effort_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'effort_cur~ effortS_prv + (1|sbj_n)');
% effort_effort_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(table_all,'effortS_cur~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_all,'effortS_cur~ effortS_prv + (1|sbj_n)');
% effortS_effortS_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% %   Hmmmm, seems effortS_prv predicts current trial effort (and effortS)...

%% Plot a regression line %
% SVtmp = bhvs{s}.SV;
% 
% figure; scatter(SVtmp,ProbAccept,'.')
% 
% sigparam=sigm_fit(SVtmp,ProbAccept,1)
% 
% [p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
% PotValues=[-10:0.1:10];
% Pout=polyval(p,PotValues);
% 
% figure; plot(PotValues,Pout)
% 
% % Test behavioral modelling.
% figure;
% subplot(2,1,1);
% scatter(bhvs{s}.stake, bhvs{s}.SV);
% title('Reward')
% subplot(2,1,2);
% scatter(bhvs{s}.effort, bhvs{s}.SV);
% title('Effort');

