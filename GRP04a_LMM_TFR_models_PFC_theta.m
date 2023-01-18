%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_madA8t1_f2t40';%'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40';%'TFRmth_D1t1_zS8t0_f2t40';%
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
outlier_thresh = 4;
n_quantiles = 5;

save_fig = 1;
fig_ftype = 'png';

if contains(an_id,'_S')
    if contains(an_id,'A8t1')
        an_lim = [-0.8 0];
    else
        an_lim = [0.5 1.5];
    end
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
out_thresh_str = ['_out' num2str(outlier_thresh)];
win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');

table_name = [an_id win_str norm_bhv_str norm_nrl_str];
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
%% Full model
% lme_full = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + dec_diff_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + dec_diff_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ effortS_cur + dec_diff_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + dec_diff_cur + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');
% 
% pfc_theta_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% pfc_theta_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% pfc_theta_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true)
% pfc_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
% pfc_theta_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
% pfc_theta_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');

pfc_theta_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
pfc_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
pfc_theta_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
pfc_theta_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_cur + SV_prv + (1|sbj_n)');

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
pfc_theta_rew_prv_vs_SV_prv = compare(lme2,lme1,'NSim',1000)
% pfc_theta_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Plot PFC theta ~ previous reward as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.PFC_theta,'reward_prv','PFC_theta',lme1,pfc_theta_rew_prv.pValue(2));
xlabel('Previous Reward (z)');
ylabel('PFC theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta ~ previous reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_theta,'reward_prv','PFC_theta',...
    lme1,pfc_theta_rew_prv.pValue(2),n_quantiles);
xlabel('Previous Reward (z)');
ylabel('PFC theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'reward','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% % Plot median splits of PFC theta as function of current/previous reward
% fig_name = 'GRP_TFR_LMM_results_PFC_theta_reward_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.PFC_theta.reward_prv<median(good_tbl_prv.PFC_theta.reward_prv);
% cur_low_idx = good_tbl_prv.PFC_theta.reward_cur<median(good_tbl_prv.PFC_theta.reward_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous Reward','High Previous Reward'});
% ylabel('PFC theta (z)');
% title('Median Split of Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('PFC theta (z)');
% title('Median Split of Low Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('PFC theta (z)');
% title('Median Split of High Previous Reward');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%%  PFC theta and previous subjective value:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

% Plot PFC theta ~ SV_prv as scatter plot 
fn_plot_LMM_scatter(SBJs,good_tbl_prv.PFC_theta,'SV_prv','PFC_theta',lme1,pfc_theta_sv_prv.pValue(2));
xlabel('Previous Subjective Value (z)');
ylabel('PFC theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta ~ previous SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_theta,'SV_prv','PFC_theta',...
    lme1,pfc_theta_sv_prv.pValue(2),n_quantiles);
xlabel('Previous Subjective Value (z)');
ylabel('PFC theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'SV','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC theta as function of increase vs. decrease in SV
%           current low   high
%   prv low                 X
%   prv high         X
% fig_name = 'GRP_TFR_LMM_results_PFC_theta_SV_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.PFC_theta.SV_prv<median(good_tbl_prv.PFC_theta.SV_prv);
% cur_low_idx = good_tbl_prv.PFC_theta.SV_cur<median(good_tbl_prv.PFC_theta.SV_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
% ylabel('PFC theta (z)');
% title('Median Split of Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('PFC theta (z)');
% title('Median Split of Low Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_theta.PFC_theta(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('PFC theta (z)');
% title('Median Split of High Previous SV');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

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

% Plot PFC theta ~ absSV as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'absSV','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'dec_diff','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


