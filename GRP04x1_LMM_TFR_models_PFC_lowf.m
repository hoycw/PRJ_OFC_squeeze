%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4';
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

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/PFC_lowf/'];
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
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all_lowf.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Toss outliers
pow_vars = {'PFC_lowf'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>st.outlier_thresh;
    
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
%   PFC THETA
%  ========================================================================
%% Full model
% lme_full = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEc = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEp = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewc = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewp = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + dec_ease_cur + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffc = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffp = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');
% 
% pfc_lowf_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% pfc_lowf_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% pfc_lowf_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true)
% pfc_lowf_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
% pfc_lowf_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
% pfc_lowf_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewc = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewp = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffc = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffp = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');

pfc_lowf_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
pfc_lowf_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
pfc_lowf_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
pfc_lowf_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_sv_curprv = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ SV_cur + SV_prv + (1|sbj_n) + (1|trl_n_cur)');
pfc_lowf_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme_sv_cur = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ SV_cur + (1|sbj_n) + (1|trl_n_cur)');
pfc_lowf_effp = compare(lme_sv_cur,lme_sv_curprv,'CheckNesting',true)

lme_ez_curprv = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ dec_ease_cur + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');

%% PFC lowf and previous reward:
lme0 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme1rs = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_prv + (1 + reward_prv|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
lme3 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ reward_prv + SV_prv + (1|sbj_n)');
pfc_lowf_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_lowf_SV_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ effort_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ effortS_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_lowf_rew_prv_vs_SV_prv = compare(lme2,lme1,'NSim',1000)
% pfc_lowf_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Plot PFC lowf ~ previous reward as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.PFC_lowf,'reward_prv','PFC_lowf',lme1,pfc_lowf_rew_prv.pValue(2));
xlabel('Previous Reward (z)');
ylabel('PFC lowf (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot lowf ~ previous reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_lowf,'reward_prv','PFC_lowf',...
    lme1,pfc_lowf_rew_prv.pValue(2),n_quantiles);
xlabel('Previous Reward (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('PFC lowf (z)');
ylim([-0.6 0.6]);
yticks(-0.5:0.5:0.5);
set(gca,'FontSize',20);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC lowf ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_lowf,'reward','PFC_lowf');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% % Plot median splits of PFC lowf as function of current/previous reward
% fig_name = 'GRP_TFR_LMM_results_PFC_lowf_reward_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.PFC_lowf.reward_prv<median(good_tbl_prv.PFC_lowf.reward_prv);
% cur_low_idx = good_tbl_prv.PFC_lowf.reward_cur<median(good_tbl_prv.PFC_lowf.reward_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous Reward','High Previous Reward'});
% ylabel('PFC lowf (z)');
% title('Median Split of Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('PFC lowf (z)');
% title('Median Split of Low Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('PFC lowf (z)');
% title('Median Split of High Previous Reward');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%% PFC lowf and decision:
% lme0 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ decision_cur + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ decision_prv + (1|sbj_n)');%,'StartMethod','random');
% pfc_lowf_dec_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_lowf_dec_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% 
% lme0 = fitglme(good_tbl_prv.PFC_lowf,'decision_cur ~ 1 + (1|sbj_n)','Distribution','binomial');
% lme1 = fitglme(good_tbl_prv.PFC_lowf,'decision_cur ~ PFC_lowf + (1|sbj_n)','Distribution','binomial');
% dec_cur_pfc_lowf = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% PFC lowf and reward change and Global Reward State (GRS):
lme0 = fitlme(good_tbl_grs.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_grs.PFC_lowf,'PFC_lowf~ reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_grs.PFC_lowf,'PFC_lowf~ grs + (1|sbj_n)');
lme3 = fitlme(good_tbl_grs.PFC_lowf,'PFC_lowf~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_lowf_rew_chg = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_lowf_grs     = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
pfc_lowf_rew_prv_vs_rew_chg = compare(lme1,lme3,'NSim',1000)
pfc_lowf_rew_prv_vs_grs     = compare(lme2,lme3,'NSim',1000)
% pfc_lowf_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

lme0 = fitlme(good_tbl_all.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_all.PFC_lowf,'PFC_lowf~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
pfc_lowf_rew_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%%  PFC lowf and previous subjective value:
lme0 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_lowf_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

% Plot PFC lowf ~ SV_prv as scatter plot 
fn_plot_LMM_scatter(SBJs,good_tbl_prv.PFC_lowf,'SV_prv','PFC_lowf',lme1,pfc_lowf_sv_prv.pValue(2));
xlabel('Previous Subjective Value (z)');
ylabel('PFC lowf (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC lowf ~ previous SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_lowf,'SV_prv','PFC_lowf',...
    lme1,pfc_lowf_sv_prv.pValue(2),n_quantiles);
xlabel('Previous Subjective Value (z)');
ylabel('PFC lowf (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC lowf ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_lowf,'SV','PFC_lowf');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot PFC lowf as function of increase vs. decrease in SV
%           current low   high
%   prv low                 X
%   prv high         X
% fig_name = 'GRP_TFR_LMM_results_PFC_lowf_SV_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.PFC_lowf.SV_prv<median(good_tbl_prv.PFC_lowf.SV_prv);
% cur_low_idx = good_tbl_prv.PFC_lowf.SV_cur<median(good_tbl_prv.PFC_lowf.SV_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
% ylabel('PFC lowf (z)');
% title('Median Split of Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('PFC lowf (z)');
% title('Median Split of Low Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.PFC_lowf.PFC_lowf(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('PFC lowf (z)');
% title('Median Split of High Previous SV');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%% PFC lowf and current reward following low reward
% prv_low_idx = good_tbl_prv.PFC_lowf.reward_prv<median(good_tbl_prv.PFC_lowf.reward_prv);
% cur_low_idx = good_tbl_prv.PFC_lowf.reward_cur<median(good_tbl_prv.PFC_lowf.reward_cur);
% lme0 = fitlme(good_tbl_prv.PFC_lowf(prv_low_idx,:),'PFC_lowf~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_lowf(prv_low_idx,:),'PFC_lowf~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
% pfc_lowf_rew_comb = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% figure; hold on;
% bins = -3:0.2:3;
% histogram(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & cur_low_idx),bins,'FaceColor','b');
% histogram(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & ~cur_low_idx),bins,'FaceColor','r');
% line([mean(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & cur_low_idx)) mean(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & cur_low_idx))],...
%     ylim,'Color','b','Linewidth',3);
% line([mean(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & ~cur_low_idx)) mean(good_tbl_prv.PFC_lowf.PFC_lowf(prv_low_idx & ~cur_low_idx))],...
%     ylim,'Color','r','Linewidth',3);
% 
%% PFC Theta Salience models
% PFC lowf salience:
lme0 = fitlme(good_tbl_all.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.PFC_lowf,'PFC_lowf~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.PFC_lowf,'PFC_lowf~ dec_ease_cur + (1|sbj_n)');
pfc_lowf_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_lowf_dec_ease = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% PFC lowf previous salience:
lme0 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.PFC_lowf,'PFC_lowf~ dec_ease_prv + (1|sbj_n)');
pfc_lowf_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_lowf_dec_ease_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Plot PFC lowf ~ absSV as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_lowf,'absSV','PFC_lowf');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
fn_plot_LMM_gratton(good_tbl_prv.PFC_lowf,'dec_ease','PFC_lowf');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


