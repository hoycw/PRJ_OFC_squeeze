%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%% Analysis parameters:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S0t2ts_wl2s025_bhvz_nrlz_out4main';
lmm_vars = {'PFC_theta','reward_prv';
            'PFC_betalo','effortS_cur';
            'BG_theta','reward_prv';
            'BG_betalo','effortS_cur'};
% lmm_vars = {'PFC_theta','SV_prv';
%             'PFC_betalo','SV_cur';
%             'BG_theta','SV_prv';
%             'BG_betalo','SV_cur'};
        
%% Analysis Set Up
% Load SBJ info, stat info:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if ~isfield(st,'time_resolved') || ~st.time_resolved; error('only use time resovled stat_ids!'); end
if st.use_simon_tfr~=0; error('why use Simon TFR?'); end

% Load power band info
[theta_cf, betalo_cf, ~] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');

%% Compute Theta and Beta power
theta_pow  = cell([numel(SBJs) 2]);
betalo_pow = cell([numel(SBJs) 2]);
bhvs       = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    if st.use_simon_tfr
        [tfr, bhvs{s}] = fn_load_simon_TFR(SBJs,st.toss_same_trials,man_trl_rej_ix);
    else
        proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
        fprintf('Loading %s\n',proc_fname);
        load(proc_fname,'tfr');
        load([sbj_dir SBJs{s} '_stim_preproc.mat']);
        bhvs{s} = sbj_data.bhv;
    end
    
    % Check channel index
    if ~any(strcmp(tfr.label{1},{'STN','GPi'})); error('BG is not first channel!'); end
    if ~any(strcmp(tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    %% Compute single trial power
    cfg = [];
    cfg.avgoverfreq = 'yes';
    cfg.avgoverrpt  = 'no';
    for ch_ix = 1:2
        cfg.channel = tfr.label(ch_ix);
        if isfield(st,'win_len')
            % Average within windows
            win_lim = fn_get_win_lim_from_center(tfr.time,st.win_center,st.win_len);
            plt_time_vec = mean(tfr.time(win_lim),2)';
            cfg.avgovertime = 'yes';
            theta_pow{s,ch_ix}  = nan([size(tfr.powspctrm,1) size(win_lim,1)]);
            betalo_pow{s,ch_ix} = nan([size(tfr.powspctrm,1) size(win_lim,1)]);
            for w_ix = 1:size(win_lim,1)
                cfg.latency     = tfr.time(win_lim(w_ix,:));
                % Theta
                cfg.frequency   = squeeze(theta_lim(s,ch_ix,:))';
                pow = ft_selectdata(cfg, tfr);
                theta_pow{s,ch_ix}(:,w_ix) = pow.powspctrm;
                % Beta Low
                cfg.frequency = squeeze(betalo_lim(s,ch_ix,:))';
                pow = ft_selectdata(cfg, tfr);
                betalo_pow{s,ch_ix}(:,w_ix) = pow.powspctrm;
            end
        else
            % Use entire time series
            cfg.avgovertime = 'no';
            cfg.latency     = st.stat_lim;
            % Theta
            cfg.frequency   = squeeze(theta_lim(s,ch_ix,:))';
            pow = ft_selectdata(cfg, tfr);
            plt_time_vec = pow.time;
            theta_pow{s,ch_ix} = squeeze(pow.powspctrm);
            % Beta Low
            cfg.frequency = squeeze(betalo_lim(s,ch_ix,:))';
            pow = ft_selectdata(cfg, tfr);
            betalo_pow{s,ch_ix} = squeeze(pow.powspctrm);
        end
    end
end

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

%% Find outliers tossed from main time-averaged analysis
% Load group model table
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' st.outlier_stat_id '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all_avg = readtable(table_all_fname);

% Identify outliers
pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
out_idx_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{f}) = abs(table_all_avg.(pow_vars{f}))>st.outlier_thresh;
    fprintf(2,'Bad trials in table_all for %s: %d\n',pow_vars{f},sum(out_idx_all.(pow_vars{f})));
end

%% Run model per time point, tossing outliers
lmm_ts      = cell([size(lmm_vars,1) length(plt_time_vec)]);
lmm_stat_ts = cell([size(lmm_vars,1) length(plt_time_vec)]);
for m_ix = 1:size(lmm_vars,1)
    fprintf('Running %d/%d LMMs: (total = %d)\n\t',m_ix,size(lmm_vars,1),length(plt_time_vec));
    for t_ix = 1:length(plt_time_vec)
        %% Create time-specific neural table
        BG_theta   = [];
        BG_betalo  = [];
        PFC_theta  = [];
        PFC_betalo = [];
        for s = 1:length(SBJs)
            BG_theta   = [BG_theta; fn_normalize_predictor(theta_pow{s,1}(:,t_ix),st.norm_nrl_pred)];
            BG_betalo  = [BG_betalo; fn_normalize_predictor(betalo_pow{s,1}(:,t_ix),st.norm_nrl_pred)];
            PFC_theta  = [PFC_theta; fn_normalize_predictor(theta_pow{s,2}(:,t_ix),st.norm_nrl_pred)];
            PFC_betalo = [PFC_betalo; fn_normalize_predictor(betalo_pow{s,2}(:,t_ix),st.norm_nrl_pred)];
        end
        
        %% Convert into table format suitable for LME modelling
        table_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, BG_theta, BG_betalo,...
            rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_ease_cur,...
            rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_ease_prv,...
            reward_chg, grs);
        
        %% Compute all desired LMMs
        % Toss outliers and NaNs
        nan_idx = isnan(table_all.(lmm_vars{m_ix,2}));
        good_tbl = table_all(~out_idx_all.(lmm_vars{m_ix,1}) & ~nan_idx,:);
        
        % Run LMM
        lme0 = fitlme(good_tbl,[lmm_vars{m_ix,1} '~ 1 + (1|sbj_n)']);
        lmm_ts{m_ix,t_ix} = fitlme(good_tbl,[lmm_vars{m_ix,1} '~ 1 + ' lmm_vars{m_ix,2} ' + (1|sbj_n)']);
        lmm_stat_ts{m_ix,t_ix} = compare(lme0,lmm_ts{m_ix,t_ix},'CheckNesting',true);%,'NSim',1000)
        
        if mod(t_ix,10)==0; fprintf('%.02f..',plt_time_vec(t_ix)); end
    end
    fprintf('\n');
end
fprintf('\n');

%% Save LMEs
if any(contains(lmm_vars(:,2),'SV'))
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMM_timeresolved_SV.mat'];
else
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMM_timeresolved.mat'];
end
fprintf('Saving %s\n',lmm_fname);
save(lmm_fname,'-v7.3','lmm_vars','lmm_ts','lmm_stat_ts','plt_time_vec');
