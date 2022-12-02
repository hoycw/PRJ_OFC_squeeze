%% Mixed modelling on the Processed data %%
!!! make sure outliers are consistently tossed!!!
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

%% Load group model tables
if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
out_thresh_str = ['out' num2str(outlier_thresh)];
table_name = [an_id out_thresh_str norm_bhv_str norm_nrl_str];

table_cur_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_cur.csv'];
fprintf('\tLoading %s...\n',table_cur_fname);
table_cur = readtable(table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_prv.csv'];
fprintf('\tLoading %s...\n',table_prv_fname);
table_prv = readtable(table_prv_fname);

% combined trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_prv_fname);
table_all = readtable(table_prv_fname);

%% Toss outliers
pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
out_idx_cur = struct;
out_idx_prv = struct;
out_ix = [];
for f = 1:length(pow_vars)
    out_idx_cur.(pow_vars{f}) = abs(table_cur.(pow_vars{f}))>outlier_thresh;
    out_idx_prv.(pow_vars{f}) = abs(table_prv.(pow_vars{f}))>outlier_thresh;
    if any(out_idx_cur.(pow_vars{f}))
        out_ix = [out_ix; find(out_idx_cur.(pow_vars{f}))];
        fprintf(2,'\t%d outliers for %s:\t',sum(out_idx_cur.(pow_vars{f})),pow_vars{f});
        fprintf(2,'%.2f, ',table_cur.(pow_vars{f})(out_idx_cur.(pow_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},outlier_thresh);
    end
end
% Toss all outliers
all_outliers = unique(out_ix);
fprintf(2,'Total bad trials: %d\n',length(all_outliers));
% table_cur(outlier_ix,:)  = [];
% table_prv(outlier_ix,:) = [];

%% ========================================================================
%   PFC THETA
%  ========================================================================
%  PFC theta and previous subjective value:
lme0 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

% Plot theta as scatter plot and median split of SV_prv
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_prv_scatter';
figure('Name',fig_name); hold on;
xvals = min(table_prv.SV_prv)-x_fudge:0.01:max(table_prv.SV_prv)+x_fudge;
for s = 1:length(SBJs)
    scatter(table_prv.SV_prv(table_prv.sbj_n==s),table_prv.PFC_theta(table_prv.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(table_prv.SV_prv(table_prv.sbj_n==s),table_prv.PFC_theta(table_prv.sbj_n==s));
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
    fig_dir   = [prj_dir 'results/TFR/LMM/' table_name '/PFC_theta/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta as function of increase vs. decrease in SV
%           current low   high
%   prv low                 X
%   prv high         X
fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_splits';
figure('Name',fig_name);
prv_low_idx = table_prv.SV_prv<median(table_prv.SV_prv);
cur_low_idx = table_all.SV_cur<median(table_all.SV_cur);

subplot(1,3,1); hold on;
vdata.lo = table_prv.PFC_theta(prv_low_idx);
vdata.hi = table_prv.PFC_theta(~prv_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
ylabel('PFC theta (z)');
title('Median Split of Previous SV');
set(gca,'FontSize',16);

subplot(1,3,2); hold on;
vdata.lo = table_all.PFC_theta(prv_low_idx & cur_low_idx);
vdata.hi = table_all.PFC_theta(prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
ylabel('PFC theta (z)');
title('Median Split of Low Previous SV');
set(gca,'FontSize',16);

subplot(1,3,3); hold on;
vdata.lo = table_all.PFC_theta(~prv_low_idx & cur_low_idx);
vdata.hi = table_all.PFC_theta(~prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
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
lL_avg = mean(table_all.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
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
lme0 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
lme3 = fitlme(table_prv(~out_idx_prv.PFC_theta,:),'PFC_theta~ reward_prv + SV_prv + (1|sbj_n)');
pfc_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_theta_SV_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(table_prv,'PFC_theta~ effort_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(table_prv,'PFC_theta~ effortS_prv + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_eff_vs_effS_prv = compare(lme1,lme2)%,'NSim',1000)
% pfc_theta_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Plot theta ~ previosu reward as scatter plot
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_prv_scatter';
figure('Name',fig_name); hold on;
xvals = min(table_prv.reward_prv)-x_fudge:0.01:max(table_prv.reward_prv)+x_fudge;
for s = 1:length(SBJs)
    scatter(table_prv.reward_prv(table_prv.sbj_n==s),table_prv.PFC_theta(table_prv.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(table_prv.reward_prv(table_prv.sbj_n==s),table_prv.PFC_theta(table_prv.sbj_n==s));
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

% Plot median splits of PFC theta as function of current/previous reward
fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_splits';
figure('Name',fig_name);
prv_low_idx = table_prv.reward_prv<median(table_prv.reward_prv);
cur_low_idx = table_all.reward_cur<median(table_all.reward_cur);

subplot(1,3,1); hold on;
vdata.lo = table_prv.PFC_theta(prv_low_idx);
vdata.hi = table_prv.PFC_theta(~prv_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Previous Reward','High Previous Reward'});
ylabel('PFC theta (z)');
title('Median Split of Previous Reward');
set(gca,'FontSize',16);

subplot(1,3,2); hold on;
vdata.lo = table_all.PFC_theta(prv_low_idx & cur_low_idx);
vdata.hi = table_all.PFC_theta(prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
ylabel('PFC theta (z)');
title('Median Split of Low Previous Reward');
set(gca,'FontSize',16);

subplot(1,3,3); hold on;
vdata.lo = table_all.PFC_theta(~prv_low_idx & cur_low_idx);
vdata.hi = table_all.PFC_theta(~prv_low_idx & ~cur_low_idx);
violins = violinplot(vdata,{'lo','hi'},'ShowMean',true);%,'ViolinAlpha',0.3);
violins(1).MeanPlot.LineWidth = 5; violins(2).MeanPlot.LineWidth = 5;
line([1 2],[violins(1).MeanPlot.YData(1) violins(2).MeanPlot.YData(1)],'Color','k','LineWidth',3);
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
lL_avg = mean(table_all.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
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
% prv_low_idx = table_all.reward_prv<median(table_all.reward_prv);
% cur_low_idx = table_all.reward_cur<median(table_all.reward_cur);
% lme0 = fitlme(table_all(prv_low_idx,:),'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(table_all(prv_low_idx,:),'PFC_theta~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_rew_comb = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% figure; hold on;
% bins = -3:0.2:3;
% histogram(table_all.PFC_theta(prv_low_idx & cur_low_idx),bins,'FaceColor','b');
% histogram(table_all.PFC_theta(prv_low_idx & ~cur_low_idx),bins,'FaceColor','r');
% line([mean(table_all.PFC_theta(prv_low_idx & cur_low_idx)) mean(table_all.PFC_theta(prv_low_idx & cur_low_idx))],...
%     ylim,'Color','b','Linewidth',3);
% line([mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx)) mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx))],...
%     ylim,'Color','r','Linewidth',3);

%% PFC Theta Salience models
% PFC theta salience:
lme0 = fitlme(table_cur,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_theta~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'PFC_theta~ dec_diff_cur + (1|sbj_n)');
pfc_theta_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_diff = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% PFC theta previous salience:
lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'PFC_theta~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'PFC_theta~ dec_diff_prv + (1|sbj_n)');
pfc_theta_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_diff_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Gratton-style line plot
fig_name = 'GRP_TFR_LMM_results_PFC_theta_salience_gratton';
figure('Name',fig_name);
subplot(1,2,1); hold on;
prv_low_idx = table_prv.absSV_prv<median(table_prv.absSV_prv);
cur_low_idx = table_all.absSV_cur<median(table_all.absSV_cur);

lL_avg = mean(table_all.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
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
prv_low_idx = table_prv.dec_diff_prv<median(table_prv.dec_diff_prv);
cur_low_idx = table_all.dec_diff_cur<median(table_all.dec_diff_cur);

lL_avg = mean(table_all.PFC_theta(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.PFC_theta(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.PFC_theta(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.PFC_theta(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.PFC_theta(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.PFC_theta(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.PFC_theta(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
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

%% ========================================================================
%   PFC BETA MODELS
%  ========================================================================
fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM_PFC_betalo/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

% PFC Beta and Reward vs. Effort models:
% PFC beta low and reward:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ reward_cur + (1|sbj_n)');
pfc_betalo_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% PFC beta low and effort:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ effort_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'PFC_betalo~ effortS_cur + (1|sbj_n)');
lme3 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1|sbj_n)');
pfc_betalo_eff = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_effS = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% pfc_betalo_effS_vs_SV = compare(lme3,lme2,'NSim',1000)
%   effortS is better model than effort
%   effort and effort S are better models than SV

% Plot theta ~ previosu reward as scatter plot
x_fudge = 0.2;
scat_sz = 20;
fig_name = 'GRP_TFR_LMM_results_PFC_betalo_effortS_scatter';
figure('Name',fig_name); hold on;
xvals = min(table_cur.effortS_cur)-x_fudge:0.01:max(table_cur.effortS_cur)+x_fudge;
for s = 1:length(SBJs)
    scatter(table_cur.effortS_cur(table_cur.sbj_n==s),table_cur.PFC_betalo(table_cur.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    mdl = fitlm(table_cur.effortS_cur(table_cur.sbj_n==s),table_cur.PFC_betalo(table_cur.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel('Previous Subj. Effort (z)');
ylabel('PFC low beta (z)');
title(['LMM beta = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(pfc_betalo_effS.pValue(2),'%.03f')]);
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Gratton-style line plot
fig_name = 'GRP_TFR_LMM_results_PFC_betalo_effort_gratton';
figure('Name',fig_name); hold on;
prv_low_idx = table_prv.effortS_prv<median(table_prv.effortS_prv);
cur_low_idx = table_all.effortS_cur<median(table_all.effortS_cur);
lL_avg = mean(table_all.PFC_betalo(prv_low_idx & cur_low_idx));
lH_avg = mean(table_all.PFC_betalo(prv_low_idx & ~cur_low_idx));
hL_avg = mean(table_all.PFC_betalo(~prv_low_idx & cur_low_idx));
hH_avg = mean(table_all.PFC_betalo(~prv_low_idx & ~cur_low_idx));
lL_sem = std(table_all.PFC_betalo(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(table_all.PFC_betalo(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(table_all.PFC_betalo(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(table_all.PFC_betalo(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel('PFC low beta (z)');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Low Previous EFF','High Previous EFF'});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{'Low Current EFF','High Current EFF'},'Location','best');
title('Effects of Previous/Current Subjective Effort');
set(gca,'FontSize',16);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% PFC beta low and previous effort:
lme0 = fitlme(table_prv,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'PFC_betalo~ effort_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'PFC_betalo~ effortS_prv + (1|sbj_n)');
pfc_betalo_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_betalo_effS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

%% ========================================================================
%   BASAL GANGLIA MODELS
%  ========================================================================

% BG beta low and subjective value:
lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betalo~ rew_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_prv,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% BG modeling section from _orig:
% fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM_BG_betalo/'];
% if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
% 
% % BG beta low and subjective value:
% lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
% lme1 = fitlme(table_cur,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)');
% bg_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % BG beta low and effort:
% lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
% lme1 = fitlme(table_cur,'BG_betalo~ effort_cur + BG_roi + (1|sbj_n)');
% lme2 = fitlme(table_cur,'BG_betalo~ effortS_cur + BG_roi + (1|sbj_n)');
% bg_betalo_eff = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% bg_betalo_effS = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% 
% % Plot betalo ~ previous reward as scatter plot
% x_fudge = 0.2;
% scat_sz = 20;
% fig_name = 'GRP_TFR_LMM_results_LFP_betalo_effortS_scatter';
% figure('Name',fig_name); hold on;
% xvals = min(table_cur.effortS_cur)-x_fudge:0.01:max(table_cur.effortS_cur)+x_fudge;
% for s = 1:length(SBJs)
%     scatter(table_cur.effortS_cur(table_cur.sbj_n==s),table_cur.BG_betalo(table_cur.sbj_n==s),...
%         scat_sz,'k');%sbj_colors(s,:));
%     mdl = fitlm(table_cur.effortS_cur(table_cur.sbj_n==s),table_cur.BG_betalo(table_cur.sbj_n==s));
%     yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
%     line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
% end
% yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
% line(xvals,yvals,'Color','k','LineWidth',5);
% xlabel('Subj. Effort (z)');
% ylabel('Basal Ganglia low beta (z)');
% title(['LMM coefficient = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(bg_betalo_effS.pValue(2),'%.03f')]);
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end
% 
% %% BG previous trial models
% % BG beta low and subjective value:
% lme0 = fitlme(table_prv,'BG_betalo~ BG_roi + (1|sbj_n)');
% lme1 = fitlme(table_prv,'BG_betalo~ SV_prv + BG_roi + (1|sbj_n)');
% bg_betalo_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% 
% % BG beta low and effort:
% lme0 = fitlme(table_prv,'BG_betalo~ BG_roi + (1|sbj_n)');
% lme1 = fitlme(table_prv,'BG_betalo~ effort_prv + BG_roi + (1|sbj_n)');
% lme2 = fitlme(table_prv,'BG_betalo~ effortS_prv + BG_roi + (1|sbj_n)');
% bg_betalo_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% bg_betalo_effS_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% 
% % Plot betalo ~ previous reward as scatter plot
% x_fudge = 0.2;
% scat_sz = 20;
% fig_name = 'GRP_TFR_LMM_results_LFP_betalo_effortS_scatter';
% figure('Name',fig_name); hold on;
% xvals = min(table_prv.effortS_prv)-x_fudge:0.01:max(table_prv.effortS_prv)+x_fudge;
% for s = 1:length(SBJs)
%     scatter(table_prv.effortS_prv(table_prv.sbj_n==s),table_prv.BG_betalo(table_prv.sbj_n==s),...
%         scat_sz,'k');%sbj_colors(s,:));
%     mdl = fitlm(table_prv.effortS_prv(table_prv.sbj_n==s),table_prv.BG_betalo(table_prv.sbj_n==s));
%     yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
%     line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
% end
% yvals = lme2.Coefficients.Estimate(1) + xvals*lme2.Coefficients.Estimate(2);
% line(xvals,yvals,'Color','k','LineWidth',5);
% xlabel('Subj. Effort (z)');
% ylabel('Basal Ganglia low beta (z)');
% title(['LMM coefficient = ' num2str(lme2.Coefficients.Estimate(2)) '; p = ' num2str(bg_betalo_effS.pValue(2),'%.03f')]);
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

% % BG theta and previous subjective value:
% lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');
% lme1 = fitlme(table_prv,'BG_theta~ reward_prv + BG_roi + (1|sbj_n)');
% bg_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% % lme0 = fitlme(table_prv,'BG_theta~ 1 + (1|sbj_n)');
% % lme1 = fitlme(table_prv,'BG_theta~ reward_prv + (1|sbj_n)');
% % bg_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)


%% ========================================================================
%   ALTERNATIVE MODELS
%  ========================================================================
%   Test Initial models:
% PFC beta low and subjective value:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1|sbj_n)');
pfc_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG beta low and subjective value:
lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_prv,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% Test alternative bands
% PFC subjective value and other bands:
lme0 = fitlme(table_cur,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_theta~ SV_cur + (1|sbj_n)');
pfc_theta_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_cur,'PFC_betahi~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betahi~ SV_cur + (1|sbj_n)');
pfc_betahi_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% PFC previous SV and other bands:
lme0 = fitlme(table_prv,'PFC_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_betalo~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betalo_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_prv,'PFC_betahi~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_betahi~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betahi_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG subjective value and other bands:
lme0 = fitlme(table_cur,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_theta~ SV_cur + (1|sbj_n)');
bg_theta_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_cur,'BG_betahi~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betahi~ SV_cur + (1|sbj_n)');
bg_betahi_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% BG previous SV and other bands:
lme0 = fitlme(table_prv,'BG_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'BG_betalo~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
bg_betalo_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
lme0 = fitlme(table_prv,'BG_betahi~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'BG_betahi~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
bg_betahi_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)


