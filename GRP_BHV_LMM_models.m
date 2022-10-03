%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_zS8t0_f2t40_log';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
log_outlier_thresh = 7;

save_fig = 1;
fig_ftype = 'png';

SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
fig_dir   = [prj_dir 'results/bhv/LMM_lRT/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

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
    for i = 1:5
        bhvs{s}.SV_effort_mn(i) = mean(bhvs{s}.SV(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_effort_se(i) = std(bhvs{s}.SV(bhvs{s}.effort==efforts(i)))./sqrt(sum(bhvs{s}.effort==efforts(i)));
        bhvs{s}.SV_stake_mn(i) = mean(bhvs{s}.SV(bhvs{s}.stake==stakes(i)));
        bhvs{s}.SV_stake_se(i) = std(bhvs{s}.SV(bhvs{s}.stake==stakes(i)))./sqrt(sum(bhvs{s}.stake==stakes(i)));
    end
end

%% Load group model tables
if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
if contains(an_id,'log'); out_thresh_str = ['out' num2str(log_outlier_thresh)]; else; out_thresh_str = ''; end

table_cur_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_cur.csv'];
fprintf('\tLoading %s...\n',table_cur_fname);
table_cur = readtable(table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_prv.csv'];
fprintf('\tLoading %s...\n',table_prv_fname);
table_prv = readtable(table_prv_fname);

% combined table
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Test log(RT) vs. RT modeling
% Reward model (this result holds for effort and SV too)
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ reward_cur + (1|sbj_n)');
rt_rew = compare(lme0,lme1)%,'NSim',1000)
lme0_log = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1_log = fitlme(table_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
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

%% Model RTs with behavioral predictors
scat_sz = 40;

fig_name = 'GRP_lRT_LMM_results';
figure('Name',fig_name);
% Reward model
lme0 = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
rt_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
subplot(2,2,1); hold on;
for s = 1:length(SBJs)
    scatter(table_cur.reward_cur(table_cur.sbj_n==s),table_cur.logrt_cur(table_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(table_cur.reward_cur):0.01:max(table_cur.reward_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Reward (z)');
xlim([-2 2]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_rew.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

% Effort model
lme0 = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'logrt_cur~ effort_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'logrt_cur~ effortS_cur + (1|sbj_n)');
rt_effort_EFF = compare(lme1,lme2)%,'NSim',1000)
%   effort is better predictor than EFFs (subjective)
rt_effort = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Subjective Value model
lme0 = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme2 = fitlme(table_cur,'logrt_cur~ SV_cur + (1|sbj_n)');
rt_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
subplot(2,2,2); hold on;
for s = 1:length(SBJs)
    scatter(table_cur.SV_cur(table_cur.sbj_n==s),table_cur.logrt_cur(table_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(table_cur.SV_cur):0.01:max(table_cur.SV_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Subjective Value (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_sv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

% Salience model
lme0 = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'logrt_cur~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'logrt_cur~ dec_diff_cur + (1|sbj_n)');
rt_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_dec_diff = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
lrt_best_diff = compare(lme1,lme2,'NSim',1000)
subplot(2,2,3); hold on;
for s = 1:length(SBJs)
    scatter(table_cur.absSV_cur(table_cur.sbj_n==s),table_cur.logrt_cur(table_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(table_cur.absSV_cur):0.01:max(table_cur.absSV_cur);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('absolute SV (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_abssv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);
subplot(2,2,4); hold on;
for s = 1:length(SBJs)
    scatter(table_cur.dec_diff_cur(table_cur.sbj_n==s),table_cur.logrt_cur(table_cur.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(table_cur.dec_diff_cur):0.01:max(table_cur.dec_diff_cur);
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Decision Difficulty (z)');
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
lme0 = fitlme(table_cur,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'logrt_cur~ decision_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'logrt_cur~ pAccept_cur + (1|sbj_n)');
rt_dec = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_pAcc = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

%% Previous trial predictors
% Previous Reward model
lme0 = fitlme(table_prv,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'logrt_cur~ reward_prv + (1|sbj_n)');
rt_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Previous Effort model
lme0 = fitlme(table_prv,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'logrt_cur~ effort_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'logrt_cur~ effortS_prv + (1|sbj_n)');
rt_effort_EFF_prv = compare(lme1,lme2)%,'NSim',1000)
rt_effort_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_effortS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Previous Subjective Value model
lme0 = fitlme(table_prv,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'logrt_cur~ SV_prv + (1|sbj_n)');
rt_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Previous salience model
lme0 = fitlme(table_prv,'logrt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'logrt_cur~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'logrt_cur~ dec_diff_prv + (1|sbj_n)');
rt_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
rt_dec_diff_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

fig_name = 'GRP_lRT_LMM_results_prv_salience';
figure('Name',fig_name);
subplot(1,2,1); hold on;
for s = 1:length(SBJs)
    scatter(table_prv.absSV_prv(table_prv.sbj_n==s),table_prv.logrt_cur(table_prv.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
end
xvals = min(table_prv.absSV_prv):0.01:max(table_prv.absSV_prv);
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Previous abs(SV) (z)');
xlim([-2.5 2.5]);
ylabel('log RT (z)');
title(['coef. = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(rt_abssv_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);
subplot(1,2,2); hold on;
for s = 1:length(SBJs)
    scatter(table_prv.dec_diff_prv(table_prv.sbj_n==s),table_prv.logrt_cur(table_prv.sbj_n==s),...
        scat_sz,sbj_colors(s,:));
    
end
xvals = min(table_prv.dec_diff_prv):0.01:max(table_prv.dec_diff_prv);
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',3);
xlabel('Previous Decision Difficulty (z)');
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
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.5 0.6]);
subplot(1,2,1); hold on;
prv_low_idx = table_all.absSV_prv<median(table_all.absSV_prv);
cur_low_idx = table_all.absSV_cur<median(table_all.absSV_cur);

lL_avg = mean(table_all.logrt_cur(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.logrt_cur(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.logrt_cur(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.logrt_cur(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.logrt_cur(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.logrt_cur(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.logrt_cur(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.logrt_cur(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('log RT');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. abs(SV)','High Prev. abs(SV)'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. abs(SV)','High Curr. abs(SV)'},'Location','northwest');
title('Effects of Previous/Current abs(SV)');
set(gca,'FontSize',16);

subplot(1,2,2); hold on;
prv_low_idx = table_all.dec_diff_prv<median(table_all.dec_diff_prv);
cur_low_idx = table_all.dec_diff_cur<median(table_all.dec_diff_cur);

lL_avg = mean(table_all.logrt_cur(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.logrt_cur(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.logrt_cur(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.logrt_cur(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.logrt_cur(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.logrt_cur(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.logrt_cur(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.logrt_cur(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('log RT');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. Difficulty','High Prev. Difficulty'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. Difficulty','High Curr. Difficulty'},'Location','northwest');
title('Effects of Previous/Current Decision Difficulty');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Does previous trial subjective value predict anything useful?
% Current subjective value
lme0 = fitlme(table_all,'SV_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'SV_cur~ SV_prv + (1|sbj_n)');
sv_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Current probability of accept
lme0 = fitlme(table_all,'pAccept_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'pAccept_cur~ SV_prv + (1|sbj_n)');
pAcc_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Current decision
lme0 = fitlme(table_all,'decision_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'decision_cur~ SV_prv + (1|sbj_n)');
dec_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Current salience
lme0 = fitlme(table_all,'absSV_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'absSV_cur~ SV_prv + (1|sbj_n)');
abssv_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_all,'dec_diff_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'dec_diff_cur~ SV_prv + (1|sbj_n)');
dec_diff_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Current reward
lme0 = fitlme(table_all,'reward_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'reward_cur~ SV_prv + (1|sbj_n)');
rew_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Current Effort
lme0 = fitlme(table_all,'effort_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'effort_cur~ SV_prv + (1|sbj_n)');
effort_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_all,'effortS_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'effortS_cur~ SV_prv + (1|sbj_n)');
effortS_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Why? check autocorrelation of effortS
lme0 = fitlme(table_all,'effort_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'effort_cur~ effortS_prv + (1|sbj_n)');
effort_effort_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_all,'effortS_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_all,'effortS_cur~ effortS_prv + (1|sbj_n)');
effortS_effortS_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
%   Hmmmm, seems effortS_prv predicts current trial effort (and effortS)...

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
