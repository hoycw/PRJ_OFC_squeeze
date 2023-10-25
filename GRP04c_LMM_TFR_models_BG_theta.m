%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4';%'S5t15_bhvz_nrl0_out4';%
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

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/BG_theta/'];
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
pow_vars = {'BG_theta'};
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
%   BASAL GANGLIA THETA
%  ========================================================================
%% Full model
% lme_full = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEc = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + dec_ease_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewc = fitlme(good_tbl_prv.BG_theta,'BG_theta~ effortS_cur + dec_ease_cur + reward_prv + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + dec_ease_cur + effortS_prv + dec_ease_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffc = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noeffp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');
% 
% bg_theta_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% bg_theta_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% bg_theta_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true)
% bg_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
% bg_theta_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
% bg_theta_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');
lme_full_norewc = fitlme(good_tbl_prv.BG_theta,'BG_theta~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');
lme_full_norewp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');
lme_full_noeffc = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');
lme_full_noeffp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');% + (1|trl_n_cur)');

bg_theta_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
bg_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
bg_theta_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
bg_theta_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

lme_full_nop = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_addprv = compare(lme_full_nop,lme_full,'CheckNesting',true)
lme_rec_rp = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_addRprv = compare(lme_full_nop,lme_rec_rp,'CheckNesting',true)

%% Test previous reward interactions
lme_full_pRcRint = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)');
lme_full_pRcEint = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:reward_prv + (1|sbj_n)');
lme_full_pRpEint = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:reward_prv + (1|sbj_n)');
bg_theta_pRcRint_p = compare(lme_full,lme_full_pRcRint,'CheckNesting',true)
bg_theta_pRcEint_p = compare(lme_full,lme_full_pRcEint,'CheckNesting',true)
bg_theta_pRpEint_p = compare(lme_full,lme_full_pRpEint,'CheckNesting',true)

lme_full_Rcpint = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_Rcpint_p = compare(lme_full,lme_full_Rcpint,'CheckNesting',true)

fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.BG_theta,'reward_prv','reward_cur','BG_theta');
ylabel('BG theta (z)');
set(gca,'FontSize',18);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');% + (1|trl_n_cur)');
lme_sv_curprv = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_cur + SV_prv + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme_sv_cur = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_cur + (1|sbj_n)');% + (1|trl_n_cur)');
lme_sv_prv = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_svc = compare(lme_sv_prv,lme_sv_curprv,'CheckNesting',true)
bg_theta_svp = compare(lme_sv_cur,lme_sv_curprv,'CheckNesting',true)

lme_ez_curprv = fitlme(good_tbl_prv.BG_theta,'BG_theta~ dec_ease_cur + dec_ease_prv + (1|sbj_n)');% + (1|trl_n_cur)');

lme_all_BG = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + BG_roi + (1|sbj_n)');% + (1|trl_n_cur)');
lme_sv_curprv_BG = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_cur + SV_prv + BG_roi + (1|sbj_n)');% + (1|trl_n_cur)');
bg_theta_full_vs_BG = compare(lme_all,lme_all_BG,'CheckNesting',true)
bg_theta_SV_vs_BG = compare(lme_sv_curprv,lme_sv_curprv_BG,'CheckNesting',true)

% Plot BG theta ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.BG_theta,'SV','BG_theta');
ylabel('BG Theta (z)');
fn_plot_LMM_gratton(good_tbl_prv.BG_theta,'reward','BG_theta');
ylabel('BG Theta (z)');
fn_plot_LMM_gratton(good_tbl_prv.BG_theta,'effortS','BG_theta');
ylabel('BG Theta (z)');

%% Plot BG theta by ROI
bg_roi_idx_all = good_tbl_all.BG_theta.BG_roi;
figure; hold on;
errorbar([1 2],[mean(good_tbl_all.BG_theta.BG_theta(bg_roi_idx_all==1)) ...
    mean(good_tbl_all.BG_theta.BG_theta(bg_roi_idx_all==2))], [std(good_tbl_all.BG_theta.BG_theta(bg_roi_idx_all==1)) ...
    std(good_tbl_all.BG_theta.BG_theta(bg_roi_idx_all==2))]);
xlabel('BG ROI');
ylabel('BG theta');
xlim([0.5 2.5]);
set(gca,'FontSize',16);

% BG theta low and reward:
lme0 = fitlme(good_tbl_all.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
bg_theta_roi = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% BG theta reward
% BG theta low and reward:
lme0 = fitlme(good_tbl_all.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_theta,'BG_theta~ reward_cur + (1|sbj_n)');
bg_theta_rew = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0_roi = fitlme(good_tbl_all.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
% lme1_roi = fitlme(good_tbl_all.BG_theta,'BG_theta~ reward_cur + BG_roi + (1|sbj_n)');
% bg_theta_roi_rew = compare(lme0_roi,lme1_roi,'CheckNesting',true)%,'NSim',1000)
% bg_theta_BGroi = compare(lme00,lme0,'CheckNesting',true)

% BG theta and previous reward:
lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + (1|sbj_n)');
bg_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
bg_theta_rew_prv_vs_SV_prv = compare(lme1,lme2,'NSim',1000)
% lme0_roi = fitlme(good_tbl_prv.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
% lme1_roi = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + BG_roi + (1|sbj_n)');
% bg_theta_roi_rew_prv = compare(lme0_roi,lme1_roi,'CheckNesting',true)%,'NSim',1000)

% Plot BG theta ~ previous reward as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.BG_theta,'reward_prv','BG_theta',lme1,bg_theta_rew_prv.pValue(2));
xlabel('Previous Reward (z)');
ylabel('BG theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot BG theta ~ previous reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.BG_theta,'reward_prv','BG_theta',...
    lme1,bg_theta_rew_prv.pValue(2),n_quantiles);
xlabel('Previous Reward (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('BG theta (z)');
ylim([-0.6 0.6]);
yticks(-0.5:0.5:0.5);
set(gca,'FontSize',20);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% BG theta and effort
lme0 = fitlme(good_tbl_all.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_theta,'BG_theta~ effortS_cur + (1|sbj_n)');
bg_theta_effS = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% BG theta and SV
% BG theta and subjective value:
lme0 = fitlme(good_tbl_all.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.BG_theta,'BG_theta~ SV_cur + (1|sbj_n)');
bg_theta_sv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0_roi = fitlme(good_tbl_all.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
% lme1_roi = fitlme(good_tbl_all.BG_theta,'BG_theta~ SV_cur + BG_roi + (1|sbj_n)');
% bg_theta_roi_sv = compare(lme0_roi,lme1_roi,'CheckNesting',true)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% lme0_roi = fitlme(good_tbl_prv.BG_theta,'BG_theta~ BG_roi + (1|sbj_n)');
% lme1_roi = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
% bg_theta_roi_sv_prv = compare(lme0_roi,lme1_roi,'CheckNesting',true)%,'NSim',1000)

% Plot BG theta ~ previous SV as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.BG_theta,'SV_prv','BG_theta',lme1,bg_theta_sv_prv.pValue(2));
xlabel('Previous Subjective Value (z)');
ylabel('BG theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot BG theta ~ previous SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.BG_theta,'SV_prv','BG_theta',...
    lme1,bg_theta_sv_prv.pValue(2),n_quantiles,'continuous',1);
xlabel('Previous Subjective Value (z)');
ylabel('BG theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% BG theta previous reward vs. SV
%   This shows they're almost identical
% lme_rew_prv = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + (1|sbj_n)');
% lme_sv_prv  = fitlme(good_tbl_prv.BG_theta,'BG_theta~ SV_prv + (1|sbj_n)');
% bg_theta_rew_sv_prv = compare(lme_rew_prv,lme_sv_prv,'NSim',1000)

%% BG theta and decision:
lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ decision_cur + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ decision_prv + (1|sbj_n)');%,'StartMethod','random');
bg_theta_dec_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
bg_theta_dec_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

lme0 = fitglme(good_tbl_prv.BG_theta,'decision_cur ~ 1 + (1|sbj_n)','Distribution','binomial');
lme1 = fitglme(good_tbl_prv.BG_theta,'decision_cur ~ BG_theta + (1|sbj_n)','Distribution','binomial');
dec_cur_bg_theta = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% BG theta and reward change and Global Reward State (GRS):
lme0 = fitlme(good_tbl_grs.BG_theta,'BG_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_grs.BG_theta,'BG_theta~ reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_grs.BG_theta,'BG_theta~ grs + (1|sbj_n)');
lme3 = fitlme(good_tbl_grs.BG_theta,'BG_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
bg_theta_rew_chg = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
bg_theta_grs     = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
bg_theta_rew_prv_vs_grs     = compare(lme2,lme3,'NSim',1000)

% Plot BG theta ~ previous SV as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_grs.BG_theta,'grs','BG_theta',lme2,bg_theta_grs.pValue(2));
xlabel('Previous Subjective Value (z)');
ylabel('BG theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot BG theta ~ previous SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_grs.BG_theta,'grs','BG_theta',...
    lme2,bg_theta_grs.pValue(2),n_quantiles,'continuous',1);
xlabel('Previous Subjective Value (z)');
ylabel('BG theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

