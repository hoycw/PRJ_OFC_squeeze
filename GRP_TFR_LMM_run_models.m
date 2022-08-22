%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
an_id = 'TFRmth_S1t2_zS8t0_f2t40_log';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
log_outlier_thresh = 7;

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
if contains(an_id,'log'); out_thresh_str = ['out' num2str(log_outlier_thresh)]; else; out_thresh_str = ''; end

table_cur_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_cur.csv'];
fprintf('\tLoading %s...\n',table_cur_fname);
table_cur = readtable(table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_prv.csv'];
fprintf('\tLoading %s...\n',table_prv_fname);
table_prv = readtable(table_prv_fname);

%% PFC theta and previous subjective value:
lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_sv_prv = compare(lme0,lme1)%,'NSim',1000)
% plotResiduals(lme1);
% plotResiduals(lme1,'fitted');
% lme1.plotPartialDependence();

%% Test Initial models:
% PFC beta low and subejctive value:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1|sbj_n)');
pfc_betalo_sv = compare(lme0,lme1)%,'NSim',1000)

% BG beta low and subejctive value:
lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_prv,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1)%,'NSim',1000)

%% Test alternative bands
% PFC subjective value and other bands:
lme0 = fitlme(table_cur,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_theta~ SV_cur + (1|sbj_n)');
pfc_theta_sv = compare(lme0,lme1)%,'NSim',1000)
lme0 = fitlme(table_cur,'PFC_betahi~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betahi~ SV_cur + (1|sbj_n)');
pfc_betahi_sv = compare(lme0,lme1)%,'NSim',1000)

% PFC previous SV and other bands:
lme0 = fitlme(table_prv,'PFC_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_betalo~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betalo_sv_prv = compare(lme0,lme1)%,'NSim',1000)
lme0 = fitlme(table_prv,'PFC_betahi~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_betahi~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_betahi_sv_prv = compare(lme0,lme1)%,'NSim',1000)

% BG subjective value and other bands:
lme0 = fitlme(table_cur,'BG_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_theta~ SV_cur + (1|sbj_n)');
bg_theta_sv = compare(lme0,lme1)%,'NSim',1000)
lme0 = fitlme(table_cur,'BG_betahi~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betahi~ SV_cur + (1|sbj_n)');
bg_betahi_sv = compare(lme0,lme1)%,'NSim',1000)

% BG previous SV and other bands:
lme0 = fitlme(table_prv,'BG_betalo~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'BG_betalo~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
bg_betalo_sv_prv = compare(lme0,lme1)%,'NSim',1000)
lme0 = fitlme(table_prv,'BG_betahi~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'BG_betahi~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
bg_betahi_sv_prv = compare(lme0,lme1)%,'NSim',1000)

%% Reward vs. Effort models:
% PFC beta low and reward:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ reward_cur + (1|sbj_n)');
pfc_betalo_rew = compare(lme0,lme1)%,'NSim',1000)

% PFC theta and previous reward:
lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
lme2 = fitlme(table_prv,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_rew_prv = compare(lme0,lme1)%,'NSim',1000)
lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_theta~ effort_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(table_prv,'PFC_theta~ effortS_prv + (1|sbj_n)');%,'StartMethod','random');
% pfc_theta_eff_vs_effS_prv = compare(lme1,lme2)%,'NSim',1000)
pfc_theta_eff_prv = compare(lme0,lme1)%,'NSim',1000)

% BG beta low and subejctive value:
lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_cur,'BG_betalo~ rew_cur + BG_roi + (1|sbj_n)');
bg_betalo_sv = compare(lme0,lme1)%,'NSim',1000)

% BG theta and previous subjective value:
lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');
lme1 = fitlme(table_prv,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');
bg_theta_sv_prv = compare(lme0,lme1)%,'NSim',1000)

%% Test salience models
% PFC theta salience:
lme0 = fitlme(table_cur,'PFC_theta~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_theta~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'PFC_theta~ dec_diff_cur + (1|sbj_n)');
pfc_theta_abssv = compare(lme0,lme1)%,'NSim',1000)
pfc_theta_dec_diff = compare(lme0,lme1)%,'NSim',1000)
