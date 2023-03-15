%% Run Circular-Linear regression on single-trial phase at each time-frequency point
% across all SBJ, separately for each regressor to predict sine and cosine
%   Uses CircularRegression from FMAtoolbox (Copyright (C) 2012 by Michaël Zugaro)
%   No random intercept for SBJ because circular phase data ranges 0 to 2*pi,
%       and the aim is to see consistent phase across SBJs
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
addpath('/Users/colinhoy/Code/Apps/FMAtoolbox/');
close all
clear all

%% Analysis parameters:
an_id = 'TFRmth_S03t2_f2t30_fourier'; stat_id = 'S5t15_bhvz_nrlphs';
reg_vars = {'PFC_theta','reward_cur';       % follow up on PFC-BG PLV ~ reward_cur finding
            'PFC_theta','reward_prv';
            'PFC_betalo','effortS_cur';
            'BG_theta','reward_cur';        % follow up on PFC-BG PLV ~ reward_cur finding
            'BG_theta','reward_prv';
            'BG_betalo','effortS_cur'};
        
%% Analysis Set Up
% Load SBJ info, stat info:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
eval(['run ' prj_dir 'scripts/PFC05_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if st.use_simon_tfr~=0; error('why use Simon TFR?'); end

% Handle single subject indexing
if length(SBJs)~=1; error('only run this for individual SBJ, probably PFC05!'); end
all_SBJs = {'PFC03','PFC04','PFC05','PFC01'};
sbj_ix = find(strcmp(all_SBJs,SBJs));

% Load power band info
[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(all_SBJs,an_id);
theta_cf = theta_cf(:,sbj_ix); betalo_cf = betalo_cf(:,sbj_ix); betahi_cf = betahi_cf(:,sbj_ix);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
beta_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');

%% Load behavioral and phase data
ang       = cell(size(SBJs));
theta_ang = cell([length(SBJs) 2]);
beta_ang  = cell([length(SBJs) 2]);
bhvs      = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load behavior and TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    bhvs{s} = sbj_data.bhv;
    
    % Check channel index
    if ~any(strcmp(tfr.label{1},{'STN','GPi'})); error('BG is not first channel!'); end
    if ~any(strcmp(tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    % Initialize data
    if s==1
        time_vec = tfr.time;
        freq_vec = tfr.freq;
    end
    
    %% Angle Extraction
    ang{s} = angle(tfr.fourierspctrm);
    
    cfg_ang = [];
    cfg_ang.avgoverfreq = 'yes';
    cfg_ang.avgovertime = 'yes';
    cfg_ang.latency     = st.stat_lim;
    for ch_ix = 1:2
        % Compute mean phase angle in T-F window
        cfg_ang.channel = tfr.label(ch_ix);
        cfg_ang.frequency = squeeze(theta_lim(s,ch_ix,:))';
        complex = ft_selectdata(cfg_ang,tfr);
        theta_ang{s,ch_ix} = squeeze(angle(complex.fourierspctrm));
        
        cfg_ang.frequency = squeeze(beta_lim(s,ch_ix,:))';
        complex = ft_selectdata(cfg_ang,tfr);
        beta_ang{s,ch_ix} = squeeze(angle(complex.fourierspctrm));
    end
end

% Concatenate angles across subjects
ang_grp = vertcat(ang{:});

%% Concatenate behavioral variables (cur = current trial, prv = previous trial)
sbj_n      = [];
PFC_roi    = [];
BG_roi     = [];
trl_n_cur  = [];

rt_cur       = [];
logrt_cur    = [];
reward_cur   = [];
effort_cur   = [];
effortS_cur  = [];
decision_cur = [];
SV_cur       = [];
absSV_cur    = [];
pAccept_cur  = [];
dec_ease_cur = [];

trl_n_prv    = [];
rt_prv       = [];
logrt_prv    = [];
reward_prv   = [];
effort_prv   = [];
effortS_prv  = [];
decision_prv = [];
SV_prv       = [];
absSV_prv    = [];
pAccept_prv  = [];
dec_ease_prv = [];

% Possible reward context predictors
reward_chg   = [];
grs          = [];

for s = 1:length(SBJs)
    % Concatenate SBJ, beta, theta values
    trl_n = size(bhvs{s}.stake,1);
    trl_n_cur = [trl_n_cur; bhvs{s}.trl];
    sbj_n = [sbj_n; num2str(ones(trl_n,1).*s)];
    if strcmp(sbj_pfc_roi{s},'OFC'); pfc_roi_ix = 1; else; pfc_roi_ix = 2; end
    PFC_roi = [PFC_roi; num2str(ones(trl_n,1).*pfc_roi_ix)];
    if strcmp(sbj_bg_roi{s},'STN'); bg_roi_ix = 1; else; bg_roi_ix = 2; end
    BG_roi = [BG_roi; num2str(ones(trl_n,1).*bg_roi_ix)];
    
    % Add behavioral variables
    rt_cur       = [rt_cur; fn_normalize_predictor(bhvs{s}.rt,st.norm_bhv_pred)];
    logrt_cur    = [logrt_cur; fn_normalize_predictor(log(bhvs{s}.rt),st.norm_bhv_pred)];
    reward_cur   = [reward_cur; fn_normalize_predictor(bhvs{s}.stake,st.norm_bhv_pred)];
    effort_cur   = [effort_cur; fn_normalize_predictor(bhvs{s}.effort,st.norm_bhv_pred)];
    effortS_cur  = [effortS_cur; fn_normalize_predictor(bhvs{s}.EFFs,st.norm_bhv_pred)];
    decision_cur = [decision_cur; fn_normalize_predictor(bhvs{s}.decision,'none')];
    SV_cur       = [SV_cur; fn_normalize_predictor(bhvs{s}.SV,st.norm_bhv_pred)];
    absSV_cur    = [absSV_cur; fn_normalize_predictor(bhvs{s}.absSV,st.norm_bhv_pred)];
    pAccept_cur  = [pAccept_cur; fn_normalize_predictor(bhvs{s}.p_accept,st.norm_bhv_pred)];
    dec_ease_cur = [dec_ease_cur; fn_normalize_predictor(bhvs{s}.dec_ease,st.norm_bhv_pred)]; % abs(p_accept - 0.5)
    
    % Add previous trial variables
    trl_n_prv    = [trl_n_prv; bhvs{s}.trl_prv];
    rt_prv       = [rt_prv; fn_normalize_predictor(bhvs{s}.rt_prv,st.norm_bhv_pred)];
    logrt_prv    = [logrt_prv; fn_normalize_predictor(log(bhvs{s}.rt_prv),st.norm_bhv_pred)];
    reward_prv   = [reward_prv; fn_normalize_predictor(bhvs{s}.stake_prv,st.norm_bhv_pred)];
    effort_prv   = [effort_prv; fn_normalize_predictor(bhvs{s}.effort_prv,st.norm_bhv_pred)];
    effortS_prv  = [effortS_prv; fn_normalize_predictor(bhvs{s}.EFFs_prv,st.norm_bhv_pred)];
    decision_prv = [decision_prv; fn_normalize_predictor(bhvs{s}.decision_prv,'none')];
    SV_prv       = [SV_prv; fn_normalize_predictor(bhvs{s}.SV_prv,st.norm_bhv_pred)];
    absSV_prv    = [absSV_prv; fn_normalize_predictor(bhvs{s}.absSV_prv,st.norm_bhv_pred)];
    pAccept_prv  = [pAccept_prv; fn_normalize_predictor(bhvs{s}.p_accept_prv,st.norm_bhv_pred)];
    dec_ease_prv = [dec_ease_prv; fn_normalize_predictor(bhvs{s}.dec_ease_prv,st.norm_bhv_pred)];
    
    % Reward context predictors
    reward_chg   = [reward_chg; fn_normalize_predictor(bhvs{s}.stake-bhvs{s}.stake_prv,st.norm_bhv_pred)];
    grs          = [grs; fn_normalize_predictor(bhvs{s}.grs,st.norm_bhv_pred)];
end

%% Convert into table format suitable for LME modelling
table_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, ...
    rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_ease_cur,...
    rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_ease_prv,...
    reward_chg, grs);

%% Find outliers tossed from main time-averaged analysis
% % Load group model table
% table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' st.outlier_stat_id '_full_table_all.csv'];
% fprintf('\tLoading %s...\n',table_all_fname);
% table_all_avg = readtable(table_all_fname);
% 
% % Identify outliers
% pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
% out_idx_all = struct;
% for f = 1:length(pow_vars)
%     % Identify outliers
%     out_idx_all.(pow_vars{f}) = abs(table_all_avg.(pow_vars{f}))>st.outlier_thresh;
%     fprintf(2,'Bad trials in table_all for %s: %d\n',pow_vars{f},sum(out_idx_all.(pow_vars{f})));
% end

%% Run group-level circular-linear regression at each time-frequency phase value
% Stats for each regressor and TFR point
betas     = nan([size(reg_vars,1) length(freq_vec) length(time_vec)]);
r2s       = nan([size(reg_vars,1) length(freq_vec) length(time_vec)]);
pvals     = nan([size(reg_vars,1) length(freq_vec) length(time_vec)]);

% Phase-Model Regression
%   CircularRegression can't do multiple regression, so one at a time
for reg_ix = 1:size(reg_vars,1)
    % Select channel and good trials
    if contains(reg_vars{reg_ix,1},'PFC'); ch_ix = 2;
    elseif contains(reg_vars{reg_ix,1},'BG'); ch_ix = 1;
    else error('must be PFC or BG'); end
    % Toss outliers and NaNs
    nan_idx = isnan(table_all.(reg_vars{reg_ix,2}));
    good_tbl = table_all(~nan_idx,:);%~out_idx_all.(reg_vars{m_ix,1}) & 
    
    fprintf('%s ~ %s (%d/%d) freq: ',reg_vars{reg_ix,1},reg_vars{reg_ix,2},reg_ix,size(reg_vars,1));
    for f_ix = 1:length(freq_vec)
        fprintf('%.3f..',freq_vec(f_ix));
        for t_ix = 1:length(time_vec)
            % function: [beta,R2,p] = CircularRegression(x,angles,p,varargin)
            %   beta = [slope, intercept]
            [slope_int, r2s(reg_ix,f_ix,t_ix), pvals(reg_ix,f_ix,t_ix)] = ...
                CircularRegression(good_tbl.(reg_vars{reg_ix,2}), ang_grp(~nan_idx,ch_ix,f_ix,t_ix));
            % Take the slope
            betas(reg_ix,f_ix,t_ix) = slope_int(1);
        end
    end
    fprintf('\n');
end

% Correct for Multiple Comparisons (regressors, times, frequencies)
[~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
fprintf('\t\t Group stats complete:');

%% Save group-level Results
grp_out_fname = [prj_dir 'data/' SBJs{1} '/' SBJs{1} '_CLreg_' an_id '_' stat_id '.mat'];
fprintf('\tSaving %s...\n',grp_out_fname);
save(grp_out_fname,'-v7.3','betas','r2s','qvals','pvals','reg_vars','time_vec','freq_vec');
