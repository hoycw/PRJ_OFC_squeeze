%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'; stat_id = 'S5t15_bhvz_nrl0_out3';%_rt21';%'S5t15_bhvz_nrl0_out4';%
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4';
% Pre-decision:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'Dn5t0_bhvz_nrlz_out4';
% Post-decision/feedback:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'D0t1_bhvz_nrlz_out4';% stat_id = 'D0t5_bhvz_nrlz_out4';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40_osr'; stat_id = 'D0t1_bhvz_nrl0_out3';

n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/PFC_theta/'];
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
pow_vars = {'PFC_theta'};
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
%   PFC THETA
%  ========================================================================
%% Full Reward-Effort model
lme_full = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_norewp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n)');
lme_full_noeffc = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_noeffp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');

pfc_theta_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
pfc_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
pfc_theta_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
pfc_theta_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

% Statistically better to add previous trial predictors?
lme_full_nop = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + (1|sbj_n)');
pfc_theta_addprv = compare(lme_full_nop,lme_full,'CheckNesting',true)
% Maybe just previous reward?
lme_rec_rp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + (1|sbj_n)');
pfc_theta_addRprv = compare(lme_full_nop,lme_rec_rp,'CheckNesting',true)

% Plot theta ~ previous reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.PFC_theta,'reward_prv','PFC_theta',...
    lme_full,pfc_theta_rewp.pValue(2),n_quantiles);
xlabel('Previous Reward (z)');
xlim([-1.8 1.8]);
xticks(-1.5:0.5:1.5);
ylabel('PFC theta (z)');
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

%% Test previous reward interactions
lme_full = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_full_pRcRint = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + (1|sbj_n)');
lme_full_pRcEint = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_cur:reward_prv + (1|sbj_n)');
lme_full_pRpEint = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + effortS_prv:reward_prv + (1|sbj_n)');
pfc_theta_pRcRint_p = compare(lme_full,lme_full_pRcRint,'CheckNesting',true)
pfc_theta_pRcEint_p = compare(lme_full,lme_full_pRcEint,'CheckNesting',true)
pfc_theta_pRpEint_p = compare(lme_full,lme_full_pRpEint,'CheckNesting',true)

lme_full_pRcRpEint = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + reward_cur:reward_prv + effortS_prv:reward_prv + (1|sbj_n)');
pfc_theta_pRcRpEint_p = compare(lme_full_pRpEint,lme_full_pRcRpEint,'CheckNesting',true)
pfc_theta_pRcRpEint_p2 = compare(lme_full_pRcRint,lme_full_pRcRpEint,'CheckNesting',true)

% [p,F,df1,df2] = anova(lme_full,'DFMethod','Satterthwaite')


% Plot partial dependence
fn_plot_LMM_interaction_partial_dependence(lme_full_pRcRint,'reward_prv','reward_cur');
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_interaction_partial_dependence(lme_full_pRpEint,'reward_prv','effortS_prv');
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_interaction_partial_dependence(lme_full_pRcEint,'reward_prv','effortS_cur');
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot partial dependence in model with both interactions
%   These look extremely similar, so I think I can trust it either way...
%   However, the interaction not included in the model no longer has a
%   gradient consistent with an interaction, which makes sense
% fn_plot_LMM_interaction_partial_dependence(lme_full_pRcRpEint,'reward_prv','reward_cur');
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '_pRcRpEint.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
% fn_plot_LMM_interaction_partial_dependence(lme_full_pRcRpEint,'reward_prv','effortS_prv');
% if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '_pRcRpEint.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot tertile line plots
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','reward_cur','PFC_theta',3,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','effortS_prv','PFC_theta',3,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','effortS_cur','PFC_theta',3,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot median splits as line plots 
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','reward_cur','PFC_theta',2,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','effortS_prv','PFC_theta',2,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end
fn_plot_LMM_quantile_line_interaction(good_tbl_prv.PFC_theta,'reward_prv','effortS_cur','PFC_theta',2,5);
if save_fig; fig_name = get(gcf,'Name'); fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname); saveas(gcf,fig_fname); end

% Plot median split bar plots
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.PFC_theta,'reward_prv','reward_cur','PFC_theta');
ylabel('PFC theta (z)');
set(gca,'FontSize',18);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
fn_plot_LMM_gratton_bar_sbj(good_tbl_prv.PFC_theta,'reward_prv','effortS_prv','PFC_theta');
ylabel('PFC theta (z)');
set(gca,'FontSize',18);
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compare reward + effort vs. SV
lme_all = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)');
lme_sv_curprv = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_cur + SV_prv + (1|sbj_n)');
pfc_theta_full_vs_SV = compare(lme_sv_curprv,lme_all,'NSim',1000)

lme_sv_cur = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_cur + (1|sbj_n)');
lme_sv_prv = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_prv + (1|sbj_n)');
pfc_theta_svc = compare(lme_sv_prv,lme_sv_curprv,'CheckNesting',true)
pfc_theta_svp = compare(lme_sv_cur,lme_sv_curprv,'CheckNesting',true)

lme_ez_curprv = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ dec_ease_cur + dec_ease_prv + (1|sbj_n)');
lme_ez_cur = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ dec_ease_cur + (1|sbj_n)');
lme_ez_prv = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ dec_ease_prv + (1|sbj_n)');
pfc_theta_ezc = compare(lme_ez_prv,lme_ez_curprv,'CheckNesting',true)
pfc_theta_ezp = compare(lme_ez_cur,lme_ez_curprv,'CheckNesting',true)

%% PFC theta and previous reward:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme1rs = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + (1 + reward_prv|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
lme3 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + SV_prv + (1|sbj_n)');
pfc_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_theta_SV_prv g= compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
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

%% PFC theta and decision:
% lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ decision_cur + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ decision_prv + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_dec_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% pfc_theta_dec_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% 
% lme0 = fitglme(good_tbl_prv.PFC_theta,'decision_cur ~ 1 + (1|sbj_n)','Distribution','binomial');
% lme1 = fitglme(good_tbl_prv.PFC_theta,'decision_cur ~ PFC_theta + (1|sbj_n)','Distribution','binomial');
% dec_cur_pfc_theta = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

%% PFC theta and reward change and Global Reward State (GRS):
lme0 = fitlme(good_tbl_grs.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_grs.PFC_theta,'PFC_theta~ reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_grs.PFC_theta,'PFC_theta~ grs + (1|sbj_n)');
lme3 = fitlme(good_tbl_grs.PFC_theta,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_rew_chg = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_grs     = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
pfc_theta_rew_prv_vs_rew_chg = compare(lme1,lme3,'NSim',1000)
pfc_theta_rew_prv_vs_grs     = compare(lme2,lme3,'NSim',1000)
% pfc_theta_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

lme0 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_rew_cur = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% lme_prewXchg = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv*reward_chg + (1|sbj_n)');%,'StartMethod','random');
% lme_prew_chg = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + reward_chg + (1|sbj_n)');%,'StartMethod','random');
lme_prew_prXchg = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + reward_prv:reward_chg + (1|sbj_n)');%,'StartMethod','random');
% previous reward is significant, and interaction is almost significant

%% PFC Theta Salience models
% PFC theta salience:
lme0 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.PFC_theta,'PFC_theta~ dec_ease_cur + (1|sbj_n)');
pfc_theta_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_ease = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% PFC theta previous salience:
lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ dec_ease_prv + (1|sbj_n)');
pfc_theta_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
pfc_theta_dec_ease_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Plot PFC theta ~ absSV as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'absSV','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
fn_plot_LMM_gratton(good_tbl_prv.PFC_theta,'dec_ease','PFC_theta');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


