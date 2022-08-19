%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
an_id = 'TFRmth_S1t2_zS8t0_f2t40_log';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%

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

table_cur_fname = [prj_dir 'data/GRP/GRP_' an_id norm_bhv_str norm_nrl_str '_full_table_cur.csv'];
fprintf('\tLoading %s...\n',table_cur_fname);
table_cur = readtable(table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' an_id norm_bhv_str norm_nrl_str '_full_table_prv.csv'];
fprintf('\tLoading %s...\n',table_prv_fname);
table_prv = readtable(table_prv_fname);

%% RT modeling
% Reward model
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ reward_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'logrt_cur~ reward_cur + (1|sbj_n)');
rt_rew = compare(lme0,lme1)%,'NSim',1000)

% Effort model
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ effort_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'rt_cur~ effortS_cur + (1|sbj_n)');
rt_effort_EFF = compare(lme1,lme2)%,'NSim',1000)
rt_effort = compare(lme0,lme1)%,'NSim',1000)

% Subjective Value model
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ SV_cur + (1|sbj_n)');
rt_sv = compare(lme0,lme1)%,'NSim',1000)

% Salience model
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ absSV_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'rt_cur~ dec_diff_cur + (1|sbj_n)');
rt_abssv = compare(lme0,lme1)%,'NSim',1000)
rt_dec_diff = compare(lme0,lme2)%,'NSim',1000)

% Decision model
lme0 = fitlme(table_cur,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'rt_cur~ decision_cur + (1|sbj_n)');
lme2 = fitlme(table_cur,'rt_cur~ pAccept_cur + (1|sbj_n)');
rt_dec = compare(lme0,lme1)%,'NSim',1000)
rt_pAcc = compare(lme0,lme2)%,'NSim',1000)

%% Previous trial predictors
% Previous Reward model
lme0 = fitlme(table_prv,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'rt_cur~ reward_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'logrt_cur~ reward_prv + (1|sbj_n)');
rt_rew_prv = compare(lme0,lme1)%,'NSim',1000)

% Previous Effort model
lme0 = fitlme(table_prv,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'rt_cur~ effort_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'rt_cur~ effortS_prv + (1|sbj_n)');
rt_effort_EFF_prv = compare(lme1,lme2)%,'NSim',1000)
rt_effort_prv = compare(lme0,lme1)%,'NSim',1000)
rt_effortS_prv = compare(lme0,lme2)%,'NSim',1000)

% Previous Subjective Value model
lme0 = fitlme(table_prv,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'rt_cur~ SV_prv + (1|sbj_n)');
rt_sv_prv = compare(lme0,lme1)%,'NSim',1000)

% Previous salience model
lme0 = fitlme(table_prv,'rt_cur~ 1 + (1|sbj_n)');
lme1 = fitlme(table_prv,'rt_cur~ absSV_prv + (1|sbj_n)');
lme2 = fitlme(table_prv,'rt_cur~ dec_diff_prv + (1|sbj_n)');
rt_abssv_prv = compare(lme0,lme1)%,'NSim',1000)
rt_dec_diff_prv = compare(lme0,lme2)%,'NSim',1000)




%% Plot a regression line %
SVtmp = bhvs{s}.SV;

figure; scatter(SVtmp,ProbAccept,'.')

sigparam=sigm_fit(SVtmp,ProbAccept,1)

[p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
PotValues=[-10:0.1:10];
Pout=polyval(p,PotValues);

figure; plot(PotValues,Pout)

% Test behavioral modelling.
figure;
subplot(2,1,1);
scatter(bhvs{s}.stake, bhvs{s}.SV);
title('Reward')
subplot(2,1,2);
scatter(bhvs{s}.effort, bhvs{s}.SV);
title('Effort');

