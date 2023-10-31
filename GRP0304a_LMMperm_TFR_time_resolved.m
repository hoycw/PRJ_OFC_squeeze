%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all
error_str = ['Why actually run this? the time-resolved LMMs are used to visualize ' ...
    'the coefficients over time, but not for significance. Thus, no need to run ' ...
    'these LMM permutations at each time point.'];
error(error_str);

%% Analysis parameters:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S0t2ts_wl2s025_bhvz_nrl0_out3main_bt1k';
pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
predictors = {'reward_cur','effortS_cur','reward_prv','effortS_prv'};

%% Analysis Set Up
% Load SBJ info, stat info:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if ~isfield(st,'nboots'); error('only run this for permutation bootstrap stat_ids!'); end
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
out_idx_all = struct;
for nrl_ix = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{nrl_ix}) = abs(table_all_avg.(pow_vars{nrl_ix}))>mean(table_all_avg.(pow_vars{nrl_ix}))...
                                    +(st.outlier_thresh*std(table_all_avg.(pow_vars{nrl_ix})));
    fprintf(2,'Bad trials in table_all for %s: %d\n',pow_vars{nrl_ix},sum(out_idx_all.(pow_vars{nrl_ix})));
end

%% Run model per time point, tossing outliers
tic;
fit_method = 'ML';
lmm_ts      = cell([length(pow_vars) length(plt_time_vec)]);
lmm_stat_ts = cell([length(pow_vars) length(plt_time_vec) length(predictors)]);
pred_pvals  = nan([length(plt_time_vec) length(predictors)]);
for nrl_ix = 1:length(pow_vars)
    fprintf('Running %d/%d LMMs: (total = %d), time elapsed = %.2f min\n\t',...
        nrl_ix,length(pow_vars),length(plt_time_vec),toc/60);
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
        tbl_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, BG_theta, BG_betalo,...
            rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_ease_cur,...
            rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_ease_prv,...
            reward_chg, grs);
        
        % Toss outliers and NaNs from previous table
        prv_nan_idx = isnan(tbl_all.SV_prv);
        good_tbl = tbl_all(~out_idx_all.(pow_vars{nrl_ix}) & ~prv_nan_idx,:);
        tbl_fields = good_tbl.Properties.VariableNames;
        for f = 1:length(tbl_fields)
            if any(isnan(good_tbl.(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
                error(['NaN is good_tbl.' tbl_fields{f}]);
            end
        end
                
        %% Compute all desired LMMs
        % Run LMM
        full_formula = [pow_vars{nrl_ix} '~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n)'];
        lmm_ts{nrl_ix,t_ix} = fitlme(good_tbl,full_formula,'FitMethod',fit_method);
        
        for p_ix = 1:length(predictors)
            null_formula = strrep(full_formula,[predictors{p_ix} ' + '],'');
            pred_lme  = fitlme(good_tbl,null_formula,'FitMethod',fit_method);
            lmm_stat_ts{nrl_ix,t_ix,p_ix} = compare(pred_lme,lmm_ts{nrl_ix,t_ix},'CheckNesting',true);
        end
        
        % Run permutations
        null_pred_pvals = nan([length(predictors) st.nboots]);
        par_tbl = good_tbl; par_nboots = st.nboots;
        parfor b_ix = 1:par_nboots
            [~, null_pred_pvals(:,b_ix)] = fn_run_LMM_null_permutation(...
                par_tbl,predictors,full_formula,fit_method);
        end
        
        % Compute p value
        for p_ix = 1:length(predictors)
            pred_pvals(t_ix,p_ix) = sum(null_pred_pvals(p_ix,:)<=lmm_stat_ts{nrl_ix,t_ix,p_ix}.pValue(2))/st.nboots;
        end
        
    end
    fprintf('\n');
end
fprintf('DONE! time elapsed = %.2f min\n',toc/60);

%% Save LMEs
if any(contains(predictors,'SV'))
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMMperm_timeresolved_SV.mat'];
else
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMMperm_timeresolved.mat'];
end
fprintf('Saving %s\n',lmm_fname);
save(lmm_fname,'-v7.3','pow_vars','predictors','lmm_ts','lmm_stat_ts','pred_pvals','plt_time_vec');
