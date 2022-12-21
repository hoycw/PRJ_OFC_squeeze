%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

% Analysis parameters:
%   an_id used for SBJ-specific band limits
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40';% an_id = 'TFRmth_D1t1_zS8t0_f2t40';
conn_metric = 'PLV';
evnt_id     = 'S';
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end
freq_ch = 'PFC';    % which channel's SBJ-specific frequency band should be used? 'PFC' or 'BG'

% Model parameters:
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%

%% Analysis Set Up
[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
betahi_lim = fn_compute_freq_lim(SBJs,betahi_cf,'betahi');

if strcmp(freq_ch,'BG')
    pk_frq_ch_ix = 1;
elseif strcmp(freq_ch,'PFC')
    pk_frq_ch_ix = 2;
else
    error('unknown peak frequency channel');
end

if ~strcmp(norm_bhv_pred,'none')
    norm_bhv_str = ['_bhv' norm_bhv_pred];
else
    norm_bhv_str = '';
end
if ~strcmp(norm_nrl_pred,'none')
    norm_nrl_str = ['_nrl' norm_nrl_pred];
else
    norm_nrl_str = '';
end

win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');
table_name = [conn_metric '_' an_id win_str norm_bhv_str norm_nrl_str];

%% Compute Theta and Beta power
theta_conn_sbj  = cell(size(SBJs));
betalo_conn_sbj = cell(size(SBJs));
betahi_conn_sbj = cell(size(SBJs));
bhvs        = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
%     proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
%     fprintf('Loading %s\n',proc_fname);
%     load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    bhvs{s} = sbj_data.bhv;
    
    % Check channel index
    if ~strcmp(sbj_data.ts.label{1},'LFP'); error('BG LFP is not first channel!'); end
    if ~any(strcmp(sbj_data.ts.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    % For decision-locked, re-align to button press
    if strcmp(evnt_id,'D')
        cfg_realign = [];
        cfg_realign.offset = round(-bhvs{s}.rt*filt_data.fsample);
        input = ft_redefinetrial(cfg_realign, sbj_data.ts);
    else
        input = sbj_data.ts;
    end
    
    % Bandpass filter data and extract power, angle, or complex data
    cfgpp = [];
    cfgpp.demean      = 'yes';
    cfgpp.hpfilter    = 'yes';
    cfgpp.hpfreq      = theta_lim(s,pk_frq_ch_ix,1);
    cfgpp.lpfilter    = 'yes';
    cfgpp.lpfreq      = theta_lim(s,pk_frq_ch_ix,2);
    if strcmp(conn_metric,'ampcorr')
        cfgpp.hilbert = 'abs';
    elseif strcmp(conn_metric,'PLV')
        cfgpp.hilbert = 'angle';
    elseif strcmp(conn_metric,'coh')
        cfgpp.hilbert = 'complex';
    end
    filt_data = ft_preprocessing(cfgpp, input);
    
    % Trim to analysis limits (after filtering to avoid edge effects)
    cfgs = [];
    cfgs.latency = an_lim;
    conn_data = ft_selectdata(cfgs,filt_data);
    
    % Compute theta connectivity
    theta_conn_sbj{s} = fn_single_trial_connectivity(conn_data,conn_metric);
    
    % Repeat for low beta
    cfgpp.lpfreq      = betalo_lim(s,pk_frq_ch_ix,1);
    cfgpp.hpfreq      = betalo_lim(s,pk_frq_ch_ix,2);
    filt_data = ft_preprocessing(cfgpp, sbj_data.ts);
    conn_data = ft_selectdata(cfgs, filt_data);
    betalo_conn_sbj{s} = fn_single_trial_connectivity(conn_data,conn_metric);
    
    % Repeat for high beta
    cfgpp.lpfreq      = betahi_lim(s,pk_frq_ch_ix,1);
    cfgpp.hpfreq      = betahi_lim(s,pk_frq_ch_ix,2);
    filt_data = ft_preprocessing(cfgpp, sbj_data.ts);
    conn_data = ft_selectdata(cfgs, filt_data);
    betahi_conn_sbj{s} = fn_single_trial_connectivity(conn_data,conn_metric);
end

%% Concatenate variables (cur = current trial, prv = previous trial)
sbj_n      = [];
% PFC_roi    = [];
% BG_roi     = [];
trl_n_cur  = [];
theta_conn  = [];
betalo_conn = [];
betahi_conn = [];

rt_cur       = [];
logrt_cur    = [];
reward_cur   = [];
effort_cur   = [];
effortS_cur  = [];
decision_cur = [];
SV_cur       = [];
absSV_cur    = [];
pAccept_cur  = [];
dec_diff_cur = [];

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
dec_diff_prv = [];

for s = 1:length(SBJs)
    % Concatenate SBJ, beta, theta values
    trl_n = size(bhvs{s}.stake,1);
    trl_n_cur = [trl_n_cur; bhvs{s}.trl];
    sbj_n = [sbj_n; num2str(ones(trl_n,1).*s)];
%     if strcmp(sbj_pfc_roi{s},'OFC'); pfc_roi_ix = 1; else; pfc_roi_ix = 2; end
%     PFC_roi = [PFC_roi; num2str(ones(trl_n,1).*pfc_roi_ix)];
%     if strcmp(sbj_bg_roi{s},'STN'); bg_roi_ix = 1; else; bg_roi_ix = 2; end
%     BG_roi = [BG_roi; num2str(ones(trl_n,1).*bg_roi_ix)];
    
    theta_conn  = [theta_conn; fn_normalize_predictor(theta_conn_sbj{s},norm_nrl_pred)];
    betalo_conn = [betalo_conn; fn_normalize_predictor(betalo_conn_sbj{s},norm_nrl_pred)];
    betahi_conn = [betahi_conn; fn_normalize_predictor(betahi_conn_sbj{s},norm_nrl_pred)];
    
    % Add behavioral variables
    rt_cur       = [rt_cur; fn_normalize_predictor(bhvs{s}.rt,norm_bhv_pred)];
    logrt_cur    = [logrt_cur; fn_normalize_predictor(log(bhvs{s}.rt),norm_bhv_pred)];
    reward_cur   = [reward_cur; fn_normalize_predictor(bhvs{s}.stake,norm_bhv_pred)];
    effort_cur   = [effort_cur; fn_normalize_predictor(bhvs{s}.effort,norm_bhv_pred)];
    effortS_cur  = [effortS_cur; fn_normalize_predictor(bhvs{s}.EFFs,norm_bhv_pred)];
    decision_cur = [decision_cur; fn_normalize_predictor(bhvs{s}.decision,norm_bhv_pred)];
    SV_cur       = [SV_cur; fn_normalize_predictor(bhvs{s}.SV,norm_bhv_pred)];
    absSV_cur    = [absSV_cur; fn_normalize_predictor(bhvs{s}.absSV,norm_bhv_pred)];
    pAccept_cur  = [pAccept_cur; fn_normalize_predictor(bhvs{s}.p_accept,norm_bhv_pred)];
    dec_diff_cur = [dec_diff_cur; fn_normalize_predictor(bhvs{s}.dec_diff,norm_bhv_pred)]; % abs(p_accept - 0.5)
    
    % Add previous trial variables
    trl_n_prv    = [trl_n_prv; bhvs{s}.trl_prv];
    rt_prv       = [rt_prv; fn_normalize_predictor(bhvs{s}.rt_prv,norm_bhv_pred)];
    logrt_prv    = [logrt_prv; fn_normalize_predictor(log(bhvs{s}.rt_prv),norm_bhv_pred)];
    reward_prv   = [reward_prv; fn_normalize_predictor(bhvs{s}.stake_prv,norm_bhv_pred)];
    effort_prv   = [effort_prv; fn_normalize_predictor(bhvs{s}.effort_prv,norm_bhv_pred)];
    effortS_prv  = [effortS_prv; fn_normalize_predictor(bhvs{s}.EFFs_prv,norm_bhv_pred)];
    decision_prv = [decision_prv; fn_normalize_predictor(bhvs{s}.decision_prv,norm_bhv_pred)];
    SV_prv       = [SV_prv; fn_normalize_predictor(bhvs{s}.SV_prv,norm_bhv_pred)];
    absSV_prv    = [absSV_prv; fn_normalize_predictor(bhvs{s}.absSV_prv,norm_bhv_pred)];
    pAccept_prv  = [pAccept_prv; fn_normalize_predictor(bhvs{s}.p_accept_prv,norm_bhv_pred)];
    dec_diff_prv = [dec_diff_prv; fn_normalize_predictor(bhvs{s}.dec_diff_prv,norm_bhv_pred)];
end

%% Convert into table format suitable for LME modelling
table_all  = table(trl_n_cur, sbj_n, theta_conn, betalo_conn, betahi_conn,...
                 rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_diff_cur,...
                 rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_diff_prv);

%% Write table for R
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tSaving %s...\n',table_all_fname);
writetable(table_all,table_all_fname);

