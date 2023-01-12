%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40';% an_id = 'TFRmth_D1t1_zS8t0_f2t40';
conn_metric = 'ampcorr';
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end
freq_ch = 'PFC';    % which channel's SBJ-specific frequency band should be used? 'PFC' or 'BG'

norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'fishz';%'none';%
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

win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');
table_name = [conn_metric '_' an_id win_str norm_bhv_str norm_nrl_str];
fig_dir   = [prj_dir 'results/TFR/LMM/' table_name out_thresh_str '/theta_conn/'];
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
conn_vars = {'theta_conn','betalo_conn','betahi_conn'};
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(conn_vars)
    % Identify outliers
    out_idx_all.(conn_vars{f}) = abs(table_all.(conn_vars{f}))>outlier_thresh;
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(conn_vars{f}) = table_all(~out_idx_all.(conn_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',conn_vars{f});
    if any(out_idx_all.(conn_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(conn_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(conn_vars{f})),conn_vars{f});
        fprintf(2,'%.2f, ',table_all.(conn_vars{f})(out_idx_all.(conn_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',conn_vars{f},outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        conn_vars{f},size(good_tbl_all.(conn_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Create current and previous trial tables
% Toss NaNs from previous table
good_tbl_prv = good_tbl_all;
for p = 1:length(conn_vars)
    prv_nan_idx = isnan(good_tbl_prv.(conn_vars{p}).SV_prv);
    good_tbl_prv.(conn_vars{p})(prv_nan_idx,:) = [];
    prv_fields = good_tbl_prv.(conn_vars{p}).Properties.VariableNames;
    for f = 1:length(prv_fields)
        if any(isnan(good_tbl_prv.(conn_vars{p}).(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
    end
end

%% ========================================================================
%   PFC THETA
%  ========================================================================
%% Full model
% lme_full = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + dec_diff_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEc = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_noDEp = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + dec_diff_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewc = fitlme(good_tbl_prv.theta_conn,'theta_conn~ effortS_cur + dec_diff_cur + reward_prv + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% lme_full_norewp = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + dec_diff_cur + effortS_prv + dec_diff_prv + (1|sbj_n) + (1|trl_n_cur)');
% 
% tconn_dec = compare(lme_full_noDEc,lme_full,'CheckNesting',true)
% tconn_dep = compare(lme_full_noDEp,lme_full,'CheckNesting',true)
% tconn_rew = compare(lme_full_norewc,lme_full,'CheckNesting',true)
% tconn_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)

%% Full model no decision ease/difficulty
lme_full = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewc = fitlme(good_tbl_prv.theta_conn,'theta_conn~ effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_norewp = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffc = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
lme_full_noeffp = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_cur + effortS_cur + reward_prv + (1|sbj_n) + (1|trl_n_cur)');

tconn_rewc = compare(lme_full_norewc,lme_full,'CheckNesting',true)
tconn_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
tconn_effc = compare(lme_full_noeffc,lme_full,'CheckNesting',true)
tconn_effp = compare(lme_full_noeffp,lme_full,'CheckNesting',true)

%% theta connectivity and previous reward:
lme0 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme0bg = fitlme(good_tbl_prv.theta_conn,'theta_conn~ BG_roi + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme1bg = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_prv + PFC_roi + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
lme3 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ reward_prv + SV_prv + (1|sbj_n)');
% theta_conn_bg = compare(lme0,lme0bg,'CheckNesting',true)%,'NSim',1000)
theta_conn_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% theta_conn_SV_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)
% lme0 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ effort_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ effortS_prv + (1|sbj_n)');%,'StartMethod','random');
% theta_conn_eff_vs_effS_prv = compare(lme1,lme2)%,'NSim',1000)
% theta_conn_eff_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)

% Plot theta connectivity ~ previous reward as scatter plot
fn_plot_LMM_scatter(SBJs,good_tbl_prv.theta_conn,'reward_prv','theta_conn',lme1,theta_conn_rew_prv.pValue(2));
xlabel('Previous Reward (z)');
ylabel('theta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta ~ previous reward as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.theta_conn,'reward_prv','theta_conn',...
    lme1,theta_conn_rew_prv.pValue(2),n_quantiles);
xlabel('Previous Reward (z)');
ylabel('theta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta connectivity ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.theta_conn,'reward','theta_conn');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% % Plot median splits of theta connectivity as function of current/previous reward
% fig_name = 'GRP_TFR_LMM_results_theta_conn_reward_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.theta_conn.reward_prv<median(good_tbl_prv.theta_conn.reward_prv);
% cur_low_idx = good_tbl_prv.theta_conn.reward_cur<median(good_tbl_prv.theta_conn.reward_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(prv_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous Reward','High Previous Reward'});
% ylabel('theta connectivity (z)');
% title('Median Split of Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('theta connectivity (z)');
% title('Median Split of Low Previous Reward');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current Reward','High Current Reward'});
% ylabel('theta connectivity (z)');
% title('Median Split of High Previous Reward');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%%  theta connectivity and previous subjective value:
lme0 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
theta_conn_sv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

% Plot theta connectivity ~ SV_prv as scatter plot 
fn_plot_LMM_scatter(SBJs,good_tbl_prv.theta_conn,'SV_prv','theta_conn',lme1,theta_conn_sv_prv.pValue(2));
xlabel('Previous Subjective Value (z)');
ylabel('theta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta connectivity ~ previous SV as line plot
fn_plot_LMM_quantile_lines(SBJs,good_tbl_prv.theta_conn,'SV_prv','theta_conn',...
    lme1,theta_conn_sv_prv.pValue(2),n_quantiles);
xlabel('Previous Subjective Value (z)');
ylabel('theta connectivity (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta connectivity ~ SV_prv as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.theta_conn,'SV','theta_conn');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% Plot theta connectivity as function of increase vs. decrease in SV
%           current low   high
%   prv low                 X
%   prv high         X
% fig_name = 'GRP_TFR_LMM_results_theta_conn_SV_splits';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.7 0.5]);
% prv_low_idx = good_tbl_prv.theta_conn.SV_prv<median(good_tbl_prv.theta_conn.SV_prv);
% cur_low_idx = good_tbl_prv.theta_conn.SV_cur<median(good_tbl_prv.theta_conn.SV_cur);
% 
% subplot(1,3,1); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(prv_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Previous SV','High Previous SV'});
% ylabel('theta connectivity (z)');
% title('Median Split of Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,2); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('theta connectivity (z)');
% title('Median Split of Low Previous SV');
% set(gca,'FontSize',16);
% 
% subplot(1,3,3); hold on;
% vdata.lo = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx & cur_low_idx);
% vdata.hi = good_tbl_prv.theta_conn.theta_conn(~prv_low_idx & ~cur_low_idx);
% violins = violinplot(vdata,{'lo','hi'});%,'ViolinAlpha',0.3);
% line([1 2],[mean(vdata.lo) mean(vdata.hi)],'Color','k','LineWidth',3);
% set(gca,'XTickLabel',{'Low Current SV','High Current SV'});
% ylabel('theta connectivity (z)');
% title('Median Split of High Previous SV');
% set(gca,'FontSize',16);
% 
% if save_fig
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%% theta connectivity and current reward following low reward
% prv_low_idx = good_tbl_prv.theta_conn.reward_prv<median(good_tbl_prv.theta_conn.reward_prv);
% cur_low_idx = good_tbl_prv.theta_conn.reward_cur<median(good_tbl_prv.theta_conn.reward_cur);
% lme0 = fitlme(good_tbl_prv.theta_conn(prv_low_idx,:),'theta_conn~ 1 + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(good_tbl_prv.theta_conn(prv_low_idx,:),'theta_conn~ reward_cur + (1|sbj_n)');%,'StartMethod','random');
% theta_conn_rew_comb = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
% figure; hold on;
% bins = -3:0.2:3;
% histogram(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & cur_low_idx),bins,'FaceColor','b');
% histogram(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & ~cur_low_idx),bins,'FaceColor','r');
% line([mean(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & cur_low_idx)) mean(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & cur_low_idx))],...
%     ylim,'Color','b','Linewidth',3);
% line([mean(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & ~cur_low_idx)) mean(good_tbl_prv.theta_conn.theta_conn(prv_low_idx & ~cur_low_idx))],...
%     ylim,'Color','r','Linewidth',3);
% 
%% Theta connectivity Salience models
% theta connectivity salience:
lme0 = fitlme(good_tbl_all.theta_conn,'theta_conn~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_all.theta_conn,'theta_conn~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(good_tbl_all.theta_conn,'theta_conn~ dec_diff_cur + (1|sbj_n)');
theta_conn_abssv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
theta_conn_dec_diff = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% theta connectivity previous salience:
lme0 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ 1 + (1|sbj_n)');
lme1 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(good_tbl_prv.theta_conn,'theta_conn~ dec_diff_prv + (1|sbj_n)');
theta_conn_abssv_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
theta_conn_dec_diff_prv = compare(lme0,lme2,'CheckNesting',true)%,'NSim',1000)

% Plot theta connectivity ~ absSV as Gratton-style line plot
fn_plot_LMM_gratton(good_tbl_prv.theta_conn,'absSV','theta_conn');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end
fn_plot_LMM_gratton(good_tbl_prv.theta_conn,'dec_diff','theta_conn');
ylabel('PFC Theta (z)');
if save_fig
    fig_name = get(gcf,'Name');
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end


