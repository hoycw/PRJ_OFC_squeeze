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
an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_madA8t1_f2t40';%'TFRmth_S1t2_zS8t0_f2t40';%
% an_id = 'TFRmth_D1t1_madS8t0_f2t40';% an_id = 'TFRmth_D1t1_zS8t0_f2t40';
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%

use_simon_tfr = 0;
toss_same_trials = 1;

if contains(an_id,'_S')
    if contains(an_id,'A8t1')
        an_lim = [-0.8 0];
    else
        an_lim = [0.5 1.5];
    end
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end

%% Analysis Set Up
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

[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
betahi_lim = fn_compute_freq_lim(SBJs,betahi_cf,'betahi');

addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

win_str = ['_' num2str(an_lim(1)) 't' num2str(an_lim(2))];
win_str = strrep(strrep(win_str,'-','n'),'.','');
table_name = [an_id win_str norm_bhv_str norm_nrl_str];

%% Compute Theta and Beta power
theta_pow  = cell([numel(SBJs) 2]);
betalo_pow = cell([numel(SBJs) 2]);
betahi_pow = cell([numel(SBJs) 2]);
bhvs       = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    if use_simon_tfr
        [tfr, bhvs{s}] = fn_load_simon_TFR(SBJs,toss_same_trials,man_trl_rej_ix);
    else
        proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
        fprintf('Loading %s\n',proc_fname);
        load(proc_fname,'tfr');
        load([sbj_dir SBJs{s} '_stim_preproc.mat']);
        bhvs{s} = sbj_data.bhv;
    end
    
    % Initialize data
    time_vec = tfr.time;
    freq_vec = tfr.freq;
    
    % Check channel index
    if ~strcmp(tfr.label{1},'LFP'); error('BG LFP is not first channel!'); end
    if ~any(strcmp(tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    %% Compute single trial power
    for ch_ix = 1:2
        % Theta
        cfg = [];
        cfg.channel = tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgovertime = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.latency     = an_lim;
        cfg.frequency   = squeeze(theta_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tfr);
        theta_pow{s,ch_ix} = pow.powspctrm;
        theta_cf(s,ch_ix) = pow.freq;
        if pow.freq-mean(theta_lim(s,ch_ix,:)) > 0.4
            fprintf(2,'\tWARNING: %s theta center freq is off: aim = %.02f, actual = %.02f; recomputing with 2-6 Hz...\n',...
                SBJs{s},mean(theta_lim(s,ch_ix,:)),pow.freq);
%             cfg.frequency = [2 2+theta_bw];
%             pow = ft_selectdata(cfg, tfr);
%             thetapk_pow{s,ch_ix} = pow.powspctrm;
%             theta_cf(s,ch_ix) = pow.freq;
        end
        
        % Beta Low
        cfg.frequency = squeeze(betalo_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tfr);
        betalo_pow{s,ch_ix} = pow.powspctrm;
        if pow.freq-mean(betalo_lim(s,ch_ix,:)) > 0.4
            fprintf(2,'\tWARNING: %s beta low center freq is off: aim = %.02f, actual = %.02f\n',...
                SBJs{s},mean(betalo_lim(s,ch_ix,:)),pow.freq);
        end
        
        % Beta High
        cfg.frequency = squeeze(betahi_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tfr);
        betahi_pow{s,ch_ix} = pow.powspctrm;
        if pow.freq-mean(betahi_lim(s,ch_ix,:)) > 0.4
            fprintf(2,'\tWARNING: %s beta high center freq is off: aim = %.02f, actual = %.02f\n',...
                SBJs{s},mean(betahi_lim(s,ch_ix,:)),pow.freq);
        end
    end
end


%% Concatenate variables (cur = current trial, prv = previous trial)
sbj_n      = [];
PFC_roi    = [];
BG_roi     = [];
trl_n_cur  = [];
PFC_theta  = [];
PFC_betalo = [];
PFC_betahi = [];
BG_theta   = [];
BG_betalo  = [];
BG_betahi  = [];

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
    
    PFC_theta  = [PFC_theta; fn_normalize_predictor(theta_pow{s,2},norm_nrl_pred)];
    PFC_betalo = [PFC_betalo; fn_normalize_predictor(betalo_pow{s,2},norm_nrl_pred)];
    PFC_betahi = [PFC_betahi; fn_normalize_predictor(betahi_pow{s,2},norm_nrl_pred)];
    BG_theta   = [BG_theta; fn_normalize_predictor(theta_pow{s,1},norm_nrl_pred)];
    BG_betalo  = [BG_betalo; fn_normalize_predictor(betalo_pow{s,1},norm_nrl_pred)];
    BG_betahi  = [BG_betahi; fn_normalize_predictor(betahi_pow{s,1},norm_nrl_pred)];
    
    % Add behavioral variables
    rt_cur       = [rt_cur; fn_normalize_predictor(bhvs{s}.rt,norm_bhv_pred)];
    logrt_cur    = [logrt_cur; fn_normalize_predictor(log(bhvs{s}.rt),norm_bhv_pred)];
    reward_cur   = [reward_cur; fn_normalize_predictor(bhvs{s}.stake,norm_bhv_pred)];
    effort_cur   = [effort_cur; fn_normalize_predictor(bhvs{s}.effort,norm_bhv_pred)];
    effortS_cur  = [effortS_cur; fn_normalize_predictor(bhvs{s}.EFFs,norm_bhv_pred)];
    decision_cur = [decision_cur; fn_normalize_predictor(bhvs{s}.decision,'none')];
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
    decision_prv = [decision_prv; fn_normalize_predictor(bhvs{s}.decision_prv,'none')];
    SV_prv       = [SV_prv; fn_normalize_predictor(bhvs{s}.SV_prv,norm_bhv_pred)];
    absSV_prv    = [absSV_prv; fn_normalize_predictor(bhvs{s}.absSV_prv,norm_bhv_pred)];
    pAccept_prv  = [pAccept_prv; fn_normalize_predictor(bhvs{s}.p_accept_prv,norm_bhv_pred)];
    dec_diff_prv = [dec_diff_prv; fn_normalize_predictor(bhvs{s}.dec_diff_prv,norm_bhv_pred)];
    
    % Reward context predictors
    reward_chg   = [reward_chg; fn_normalize_predictor(bhvs{s}.stake-bhvs{s}.stake_prv,norm_bhv_pred)];
    grs          = [grs; fn_normalize_predictor(bhvs{s}.grs,norm_bhv_pred)];
end

%% Convert into table format suitable for LME modelling
table_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
                 rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_diff_cur,...
                 rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_diff_prv,...
                 reward_chg, grs);

%% Write table for R
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tSaving %s...\n',table_all_fname);
writetable(table_all,table_all_fname);

