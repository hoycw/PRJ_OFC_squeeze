%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%% Analysis parameters:
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'ampcorr_Sn8t0_bhvz_nrlfz_out4';
% Stimulus decision phase:
% an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'ampcorr_S5t15_bhvz_nrlfz_out4';
% Pre-decision:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'ampcorr_Dn5t0_bhvz_nrlfz_out4';
% Post-decision/feedback:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'ampcorr_D0t5_bhvz_nrlfz_out4';

% Phase-locking value
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'PLV_S5t15_bhvz_nrlz_out4';

% Jackknife coherence
% an_id = 'TFRmth_S03t2_f2t30_fourier'; stat_id = 'cohjk_S5t15_bhvz_nrlz';

%% Analysis Set Up
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/PFC05_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);

% Handle single subject indexing
if length(SBJs)~=1; error('only run this for individual SBJ, probably PFC05!'); end
all_SBJs = {'PFC03','PFC04','PFC05','PFC01'};
sbj_ix = find(strcmp(all_SBJs,SBJs));

% Load power band info
[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(all_SBJs,an_id);
theta_cf = theta_cf(:,sbj_ix); betalo_cf = betalo_cf(:,sbj_ix); betahi_cf = betahi_cf(:,sbj_ix);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');

if strcmp(st.freq_pk_ch,'BG')
    pk_frq_ch_ix = 1;
elseif strcmp(st.freq_pk_ch,'PFC')
    pk_frq_ch_ix = 2;
else
    error('unknown peak frequency channel');
end

%% Compute Theta and Beta power
theta_conn_sbj  = cell(size(SBJs));
betalo_conn_sbj = cell(size(SBJs));
bhvs        = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    bhvs{s} = sbj_data.bhv;
    
    % Check channel index
    if ~any(strcmp(sbj_data.ts.label{1},{'STN','GPi'})); error('BG is not first channel!'); end
    if ~any(strcmp(sbj_data.ts.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    if strcmp(st.conn_metric,'cohjk')
        proc_fname = [sbj_dir SBJs{s} '_' an_id '_' st.conn_metric '.mat'];
        fprintf('Loading %s\n',proc_fname);
        load(proc_fname,'conn');
        
        % Check channel index
        if ~any(strcmp(conn.label{1},{'STN','GPi'})); error('BG is not first channel!'); end
        if ~any(strcmp(conn.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
        
        % Find frequency range
        [~,t_lo_ix] = min(abs(conn.freq-theta_lim(s,pk_frq_ch_ix,1)));
        [~,t_hi_ix] = min(abs(conn.freq-theta_lim(s,pk_frq_ch_ix,2)));
        [~,b_lo_ix] = min(abs(conn.freq-betalo_lim(s,pk_frq_ch_ix,1)));
        [~,b_hi_ix] = min(abs(conn.freq-betalo_lim(s,pk_frq_ch_ix,2)));
        
        % Average across frequencies then time
        theta_conn_sbj{s}  = squeeze(nanmean(nanmean(conn.cohspctrm(:,1,2,t_lo_ix:t_hi_ix,:),4),5));
        betalo_conn_sbj{s} = squeeze(nanmean(nanmean(conn.cohspctrm(:,1,2,b_lo_ix:b_hi_ix,:),4),5));
    else
        % For decision-locked, re-align to button press
        if strcmp(st.stat_evnt,'D')
            cfg_realign = [];
            cfg_realign.offset = round(-bhvs{s}.rt*sbj_data.ts.fsample);
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
        if strcmp(st.conn_metric,'ampcorr')
            cfgpp.hilbert = 'abs';
        elseif strcmp(st.conn_metric,'PLV')
            cfgpp.hilbert = 'angle';
        elseif strcmp(st.conn_metric,'coh')
            cfgpp.hilbert = 'complex';
        end
        filt_data = ft_preprocessing(cfgpp, input);
        
        % Trim to analysis limits (after filtering to avoid edge effects)
        cfgs = [];
        cfgs.latency = st.stat_lim;
        conn_data = ft_selectdata(cfgs,filt_data);
        
        % Compute theta connectivity
        theta_conn_sbj{s} = fn_connectivity_single_trial(conn_data,st.conn_metric);
        
        % Repeat for low beta
        cfgpp.lpfreq      = betalo_lim(s,pk_frq_ch_ix,1);
        cfgpp.hpfreq      = betalo_lim(s,pk_frq_ch_ix,2);
        filt_data = ft_preprocessing(cfgpp, sbj_data.ts);
        conn_data = ft_selectdata(cfgs, filt_data);
        betalo_conn_sbj{s} = fn_connectivity_single_trial(conn_data,st.conn_metric);
        
%         % Repeat for high beta
%         cfgpp.lpfreq      = betahi_lim(s,pk_frq_ch_ix,1);
%         cfgpp.hpfreq      = betahi_lim(s,pk_frq_ch_ix,2);
%         filt_data = ft_preprocessing(cfgpp, sbj_data.ts);
%         conn_data = ft_selectdata(cfgs, filt_data);
%         betahi_conn_sbj{s} = fn_connectivity_single_trial(conn_data,st.conn_metric);
    end
end

%% Concatenate variables (cur = current trial, prv = previous trial)
sbj_n      = [];
PFC_roi    = [];
BG_roi     = [];
trl_n_cur  = [];
theta_conn  = [];
betalo_conn = [];

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
    
    theta_conn  = [theta_conn; fn_normalize_predictor(theta_conn_sbj{s},st.norm_nrl_pred)];
    betalo_conn = [betalo_conn; fn_normalize_predictor(betalo_conn_sbj{s},st.norm_nrl_pred)];
    
    % Add behavioral variables
    rt_cur       = [rt_cur; fn_normalize_predictor(bhvs{s}.rt,st.norm_bhv_pred)];
    logrt_cur    = [logrt_cur; fn_normalize_predictor(log(bhvs{s}.rt),st.norm_bhv_pred)];
    reward_cur   = [reward_cur; fn_normalize_predictor(bhvs{s}.stake,st.norm_bhv_pred)];
    effort_cur   = [effort_cur; fn_normalize_predictor(bhvs{s}.effort,st.norm_bhv_pred)];
    effortS_cur  = [effortS_cur; fn_normalize_predictor(bhvs{s}.EFFs,st.norm_bhv_pred)];
    decision_cur = [decision_cur; fn_normalize_predictor(bhvs{s}.decision,st.norm_bhv_pred)];
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
    decision_prv = [decision_prv; fn_normalize_predictor(bhvs{s}.decision_prv,st.norm_bhv_pred)];
    SV_prv       = [SV_prv; fn_normalize_predictor(bhvs{s}.SV_prv,st.norm_bhv_pred)];
    absSV_prv    = [absSV_prv; fn_normalize_predictor(bhvs{s}.absSV_prv,st.norm_bhv_pred)];
    pAccept_prv  = [pAccept_prv; fn_normalize_predictor(bhvs{s}.p_accept_prv,st.norm_bhv_pred)];
    dec_ease_prv = [dec_ease_prv; fn_normalize_predictor(bhvs{s}.dec_ease_prv,st.norm_bhv_pred)];
    
    % Reward context predictors
    reward_chg   = [reward_chg; fn_normalize_predictor(bhvs{s}.stake-bhvs{s}.stake_prv,st.norm_bhv_pred)];
    grs          = [grs; fn_normalize_predictor(bhvs{s}.grs,st.norm_bhv_pred)];
end

%% Convert into table format suitable for LME modelling
table_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, theta_conn, betalo_conn,...
                 rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_ease_cur,...
                 rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_ease_prv,...
                 reward_chg, grs);

%% Write table for R
table_all_fname = [prj_dir 'data/' SBJs{1} '/' SBJs{1} '_' an_id '_' stat_id '_full_table_all.csv'];
fprintf('\tSaving %s...\n',table_all_fname);
writetable(table_all,table_all_fname);

