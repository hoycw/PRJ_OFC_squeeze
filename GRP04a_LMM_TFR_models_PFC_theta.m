%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
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
fig_dir   = [prj_dir 'results/TFR/LMM/' table_name out_thresh_str '/PFC_theta/'];
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
%   PFC THETA
%  ========================================================================
%  PFC theta and previous subjective value:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

% Plot theta as scatter plot and median split of SV_prv
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_prv_scatter';
figure('Name',fig_name); hold on;
xvals = min(good_tbl_prv.PFC_theta.SV_prv)-x_fudge:0.01:max(good_tbl_prv.PFC_theta.SV_prv)+x_fudge;
for s = 1:length(SBJs)
    scatter(good_tbl_prv.PFC_theta.SV_prv(good_tbl_prv.PFC_theta.sbj_n==s),good_tbl_prv.PFC_theta.PFC_theta(good_tbl_prv.PFC_theta.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(good_tbl_prv.PFC_theta.SV_prv(good_tbl_prv.PFC_theta.sbj_n==s),good_tbl_prv.PFC_theta.PFC_theta(good_tbl_prv.PFC_theta.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel('Previous SV (z)');
ylabel('PFC theta (z)');
title(['LMM beta = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(pfc_theta_sv_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta as function of increase vs. decrease in SV
%           current low   high
%   prv low                 X
%   prv high         X
fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_splits';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
prv_low_idx = good_tbl_prv.PFC_theta.SV_prv<median(good_tbl_prv.PFC_theta.SV_prv);
cur_low_idx = good_tbl_prv.PFC_theta.SV_cur<median(good_tbl_prv.PFC_theta.SV_cur);

subplot(1,3,1); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
ylabel('PFC theta (z)');
title('Median Split of Previous SV');
set(gca,'FontSize',16);

subplot(1,3,2); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
ylabel('PFC theta (z)');
title('Median Split of Low Previous SV');
set(gca,'FontSize',16);

subplot(1,3,3); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
ylabel('PFC theta (z)');
title('Median Split of High Previous SV');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Gratton-style line plot
fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_gratton';
figure('Name',fig_name); hold on;
lL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('PFC Theta (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Current SV','High Current SV'},'Location','best');
title('Effects of Previous and Current Subjective Value');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% PFC theta and previous reward:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
lme3 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + SV_prv + (1|sbj_n)');
pfc_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_theta_SV_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ effort_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ effortS_prv + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_eff_vs_effS_prv = compare(lme1,lme2)%,'NSim',1000)
% pfc_theta_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Plot theta ~ previous reward as scatter plot
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_prv_scatter';
figure('Name',fig_name); hold on;
xvals = min(good_tbl_prv.PFC_theta.reward_prv)-x_fudge:0.01:max(good_tbl_prv.PFC_theta.reward_prv)+x_fudge;
for s = 1:length(SBJs)
    scatter(good_tbl_prv.PFC_theta.reward_prv(good_tbl_prv.PFC_theta.sbj_n==s),good_tbl_prv.PFC_theta.PFC_theta(good_tbl_prv.PFC_theta.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(good_tbl_prv.PFC_theta.reward_prv(good_tbl_prv.PFC_theta.sbj_n==s),good_tbl_prv.PFC_theta.PFC_theta(good_tbl_prv.PFC_theta.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel('Previous Reward (z)');
ylabel('PFC theta (z)');
title(['LMM beta = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(pfc_theta_rew_prv.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta ~ previous reward as line plot
for s = 1:length(SBJs)
    % Get bin edges
    qs = quantile(good_tbl_prv.PFC_theta.reward_prv(good_tbl_prv.PFC_theta.sbj_n==s),n_quantiles);
    rew_prv_edges = nan([n_quantiles-1 1]);
    for i = 1:n_quantiles-1
        rew_prv_edges(i) = mean(qs(i:i+1));
    end
    % Average within bins
    rew_prv_mn.(SBJs{s})            = nan([n_quantiles 1]);
    PFC_theta_rew_prv_mn.(SBJs{s})  = nan([n_quantiles 1]);
    PFC_theta_rew_prv_se.(SBJs{s})  = nan([n_quantiles 1]);
    for i = 1:n_quantiles
        if i==1                 % first quantile
            trl_idx = good_tbl_prv.PFC_theta.sbj_n==s & good_tbl_prv.PFC_theta.reward_prv<rew_prv_edges(i);
        elseif i==n_quantiles   % last quantile
            trl_idx = good_tbl_prv.PFC_theta.sbj_n==s & good_tbl_prv.PFC_theta.reward_prv>=rew_prv_edges(i-1);
        else                    % middle quantiles
            trl_idx = good_tbl_prv.PFC_theta.sbj_n==s & good_tbl_prv.PFC_theta.reward_prv<rew_prv_edges(i) &...
                good_tbl_prv.PFC_theta.reward_prv>=rew_prv_edges(i-1);
        end
        rew_prv_mn.(SBJs{s})(i) = mean(good_tbl_prv.PFC_theta.reward_prv(trl_idx));
        PFC_theta_rew_prv_mn.(SBJs{s})(i) = mean(good_tbl_prv.PFC_theta.PFC_theta(trl_idx));
        PFC_theta_rew_prv_se.(SBJs{s})(i) = std(good_tbl_prv.PFC_theta.PFC_theta(trl_idx))./sqrt(sum(trl_idx));
    end
end
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_prv_line';
figure('Name',fig_name); hold on;
xvals = min(good_tbl_prv.PFC_theta.reward_prv)-x_fudge:0.01:max(good_tbl_prv.PFC_theta.reward_prv)+x_fudge;
for s = 1:length(SBJs)
    errorbar(rew_prv_mn.(SBJs{s}),PFC_theta_rew_prv_mn.(SBJs{s}),PFC_theta_rew_prv_se.(SBJs{s}),...
        'Color',sbj_colors(s,:),'LineWidth',1.5);
end
yvals = lme1.Coefficients.Estimate(1) + xvals*lme1.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',4);
xlabel('Previous Reward (z)');
ylabel('PFC theta (z)');
title(['LMM beta = ' num2str(lme1.Coefficients.Estimate(2)) '; p = ' num2str(pfc_theta_rew_prv.pValue(2),'%.03f')]);
legend([SBJs 'Group'],'Location','best');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot median splits of PFC theta as function of current/previous reward
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_splits';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
prv_low_idx = good_tbl_prv.PFC_theta.reward_prv<median(good_tbl_prv.PFC_theta.reward_prv);
cur_low_idx = good_tbl_prv.PFC_theta.reward_cur<median(good_tbl_prv.PFC_theta.reward_cur);

subplot(1,3,1); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Previous Reward','High Previous Reward'});
ylabel('PFC theta (z)');
title('Median Split of Previous Reward');
set(gca,'FontSize',16);

subplot(1,3,2); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
ylabel('PFC theta (z)');
title('Median Split of Low Previous Reward');
set(gca,'FontSize',16);

subplot(1,3,3); hold on;
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
ylabel('PFC theta (z)');
title('Median Split of High Previous Reward');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Gratton-style line plot
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_gratton';
figure('Name',fig_name); hold on;
lL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx));
vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx);
vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx);
lL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('PFC Theta (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. Reward','High Prev. Reward'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. Reward','High Curr. Reward'},'Location','southeast');
title('Effects of Previous and Current Reward');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% PFC theta and current reward following low reward
% prv_low_idx = good_tbl_prv.PFC_theta.reward_prv<median(good_tbl_prv.PFC_theta.reward_prv);
% cur_low_idx = good_tbl_prv.PFC_theta.reward_cur<median(good_tbl_prv.PFC_theta.reward_cur);
% lme0 = fitlme(good_tbl_prv.PFC_theta(prv_low_idx,:),'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_theta(prv_low_idx,:),'PFC_theta~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_rew_comb = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% figure; hold on;
% bins = -3:0.2:3;
% histogram(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx),bins,'FaceColor','b');
% histogram(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx),bins,'FaceColor','r');
% line([mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx)) mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx))],...
%     ylim,'Color','b','Linewidth',3);
% line([mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx)) mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx))],...
%     ylim,'Color','r','Linewidth',3);
% 
%% PFC Theta Salience models
% PFC theta salience:
lme0 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ dec_diff_cur + (1|sbj_n)');
pfc_theta_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_diff = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% PFC theta previous salience:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ dec_diff_prv + (1|sbj_n)');
pfc_theta_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_diff_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Gratton-style line plot
fig_name = 'GRP_TFR_LMM_results_PFC_theta_salience_gratton';
figure('Name',fig_name);
subplot(1,2,1); hold on;
prv_low_idx = good_tbl_prv.PFC_theta.absSV_prv<median(good_tbl_prv.PFC_theta.absSV_prv);
cur_low_idx = good_tbl_prv.PFC_theta.absSV_cur<median(good_tbl_prv.PFC_theta.absSV_cur);

lL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('PFC Theta (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. abs(SV)','High Prev. abs(SV)'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. abs(SV)','High Curr. abs(SV)'},'Location','southwest');
title('Effects of Previous/Current abs(SV)');
set(gca,'FontSize',16);

subplot(1,2,2); hold on;
prv_low_idx = good_tbl_prv.PFC_theta.dec_diff_prv<median(good_tbl_prv.PFC_theta.dec_diff_prv);
cur_low_idx = good_tbl_prv.PFC_theta.dec_diff_cur<median(good_tbl_prv.PFC_theta.dec_diff_cur);

lL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('PFC Theta (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Prev. Difficulty','High Prev. Difficulty'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Curr. Difficulty','High Curr. Difficulty'},'Location','southwest');
title('Effects of Previous/Current Decision Difficulty');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

