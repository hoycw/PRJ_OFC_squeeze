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
an_id = 'TFRmth_S1t2_zS8t0_f2t40';%'TFRmth_S1t2_zS1t0_f2t40_log';%
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'logz';%'none';%
outlier_thresh = 4;
% an_id = 'TFRmth_D1t1_zS8t0_f2t40';
use_simon_tfr = 0;
toss_same_trials = 1;
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
    plt_lim = [-0.2 2];
    betahi_cf = ones([2 numel(SBJs)])*-1;
    if use_simon_tfr
        % Simon params: [10,17,13,12]
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf = [-1 -1 -1 -1; 10 17 13 12]; % PFC03, PFC04, PFC05, PFC01
    elseif strcmp(an_id,'TFRmth_S1t2_zS1t0_f2t40')  % equivalent to Simon
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; 10 17 13 12];
    elseif any(strcmp(an_id,{'TFRmth_S1t2_zS8t0_f2t40','TFRmth_S1t2_madS8t0_f2t40'}))
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; 11 19 12 12];
    elseif strcmp(an_id,'TFRw_S25t2_dbS25t05_fl2t40_c7')
        theta_cf = [2.5 3.5 3.5 6; 5 4 3.5 3]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; -1 -1 -1 -1];
    else
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf = [-1 -1 -1 -1; 10 17 13 12]; % PFC03, PFC04, PFC05, PFC01
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
%         theta_cf = [2.5 3 3.5 5; 5 3.5 3.5 3]; % BG (row1) then PFC (row2)
%         betalo_cf  = [17 22 14 14; 13 17 11 15];
        %     betahi_cf  = % no SBJ-specific peaks;
        %     betalo_cf = ones([2 numel(SBJs)])*-1;
    end
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
    betahi_cf = ones([2 numel(SBJs)])*-1;
    if strcmp(an_id,'TFRw_D1t1_dbS25t05_fl2t40_c7')
        theta_cf = [3.5 3 4.5 6; 3 3 4.5 3]; % BG (row1) then PFC (row2)
        betalo_cf  = [13 20 17 18; -1 -1 -1 -1];
        betahi_cf = ones([2 numel(SBJs)])*-1;
    elseif strcmp(an_id,'TFRmth_D1t1_zS1t0_f2t40')  % equivalent to Simon
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; -1 -1 -1 -1];%10 17 13 12];
    elseif contains(an_id,'TFRmth_D1t1_zS8t0_f2t40')
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; 11 19 12 12];
    else
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
        theta_cf = [3.5 3 3.5 7; 4.5 3 3.5 3]; % BG (row1) then PFC (row2)
        %     betalo_cf  = [17 22 14 14; 13 17 11 15];
        %     betahi_cf  = % no SBJ-specific peaks;
        betalo_cf = [20 17 11 24; -1 -1 -1 -1];
    end
end

if ~strcmp(norm_bhv_pred,'none'); norm_bhv_str = ['_bhv' norm_bhv_pred]; else; norm_bhv_str = ''; end
if ~strcmp(norm_nrl_pred,'none'); norm_nrl_str = ['_nrl' norm_nrl_pred]; else; norm_nrl_str = ''; end
out_thresh_str = ['_out' num2str(outlier_thresh)];

% Simon originals:
theta_bw  = 3;
betalo_bw = 4;
betahi_bw = 4;
theta_canon  = [4 7];
betalo_canon = [12 20];
betahi_canon = [21 35];
% theta_lim  = [4 7];
% beta_lim   = [13 30];    
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];

%% Analysis Set Up
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

theta_lim  = nan(length(SBJs),2,2);
betalo_lim = nan(length(SBJs),2,2);
betahi_lim = nan(length(SBJs),2,2);
for s = 1:length(SBJs)
    for ch_ix = 1:2
        if theta_cf(ch_ix,s)==-1
            theta_lim(s,ch_ix,:) = theta_canon;
        else
            theta_lim(s,ch_ix,:) = [theta_cf(ch_ix,s)-theta_bw/2 theta_cf(ch_ix,s)+theta_bw/2];
        end
        if betalo_cf(ch_ix,s)==-1
            betalo_lim(s,ch_ix,:) = betalo_canon;
        else
            betalo_lim(s,ch_ix,:)  = [betalo_cf(ch_ix,s)-betalo_bw/2 betalo_cf(ch_ix,s)+betalo_bw/2];
        end
        if betahi_cf(ch_ix,s)==-1
            betahi_lim(s,ch_ix,:) = betahi_canon;
        else
            betahi_lim(s,ch_ix,:)  = [betahi_cf(ch_ix,s)-betahi_bw/2 betahi_cf(ch_ix,s)+betahi_bw/2];
        end
        
        % Print final parameters
        if ch_ix==1; ch_lab = sbj_bg_roi{s}; else; ch_lab = sbj_pfc_roi{s}; end
%         fprintf('%s %s theta CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,theta_cf(ch_ix,s),theta_lim(s,ch_ix,1),theta_lim(s,ch_ix,2));
%         fprintf('%s %s beta lo CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betalo_cf(ch_ix,s),betalo_lim(s,ch_ix,1),betalo_lim(s,ch_ix,2));
%         fprintf('%s %s beta hi CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betahi_cf(ch_ix,s),betahi_lim(s,ch_ix,1),betahi_lim(s,ch_ix,2));
    end
end

%% Load Simon file details
if use_simon_tfr
    DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
    DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';
    
    [numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
    % SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
    if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end
    
    FileDetails = strings(2:end,:);
end

%% Compute Theta and Beta power
theta_pow  = cell([numel(SBJs) 2]);
betalo_pow = cell([numel(SBJs) 2]);
betahi_pow = cell([numel(SBJs) 2]);
bhvs       = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    if use_simon_tfr
        % Load data files
        proc_fname = strcat(DataStorage2,'/',FileDetails{s,1},'Stimulus_Locked.mat');
        fprintf('Loading %s\n',proc_fname);
        load(proc_fname);
        tfr = AllData.TFbl;
        bhvs{s} = AllData.exp;
        if toss_same_trials
            load([sbj_dir SBJs{s} '_stim_preproc.mat'],'sbj_data');
            % Combine bad behavioral and neural trials 
            %   (empty and bad key are already tossed in Simon data)
            simon_trl_idx = 1:150;
            if strcmp(SBJs{s},'PFC04')
                simon_trl_idx(72) = [];% from whrc neural variance
            elseif strcmp(SBJs{s},'PFC05')
                simon_trl_idx(76) = [];% from whrempty
            elseif strcmp(SBJs{s},'PFC01')
                simon_trl_idx(26) = [];% from whrtrl
            end
            all_bad_ix = unique([man_trl_rej_ix{s}'; sbj_data.bhv.empty_ix; sbj_data.bhv.bad_resp_ix; sbj_data.bhv.bad_rt_ix]);
            simon_bad_ix = [];
            for t = 1:length(all_bad_ix)
                simon_bad_ix = [simon_bad_ix; find(simon_trl_idx==all_bad_ix(t))];
            end
            if isempty(simon_bad_ix)
                fprintf('%s: All bad trials already tossed!\n',SBJs{s});
            else
                fprintf(2,'%s: Removing %d trials!\n',SBJs{s},length(simon_bad_ix));
                % Remove from behavior
                bhv_fields = fieldnames(bhvs{s});
                for f = 1:length(bhv_fields)
                    if length(bhvs{s}.(bhv_fields{f}))==length(simon_trl_idx) && ~contains(bhv_fields{f},'_ix')
                        bhvs{s}.(bhv_fields{f})(simon_bad_ix) = [];
                    end
                end
                % Remove from neural
                tfr.powspctrm(simon_bad_ix,:,:,:) = [];
            end
        end
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
    
    % Compute single trial power
    for ch_ix = 1:2
        % Theta
        cfg = [];
        cfg.channel = tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgovertime = 'no';
        cfg.avgoverrpt  = 'no';
        cfg.latency     = plt_lim;
        cfg.frequency   = squeeze(theta_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tfr);
        plt_time_vec = pow.time;
        theta_pow{s,ch_ix} = squeeze(pow.powspctrm);
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
        betalo_pow{s,ch_ix} = squeeze(pow.powspctrm);
        if pow.freq-mean(betalo_lim(s,ch_ix,:)) > 0.4
            fprintf(2,'\tWARNING: %s beta low center freq is off: aim = %.02f, actual = %.02f\n',...
                SBJs{s},mean(betalo_lim(s,ch_ix,:)),pow.freq);
        end
        
        % Beta High
        cfg.frequency = squeeze(betahi_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tfr);
        betahi_pow{s,ch_ix} = squeeze(pow.powspctrm);
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
    if strcmp(sbj_pfc_roi{s},'OFC'); pfc_roi_ix = 1; else; pfc_roi_ix = 2; end
    PFC_roi = [PFC_roi; num2str(ones(trl_n,1).*pfc_roi_ix)];
    if strcmp(sbj_bg_roi{s},'STN'); bg_roi_ix = 1; else; bg_roi_ix = 2; end
    BG_roi = [BG_roi; num2str(ones(trl_n,1).*bg_roi_ix)];
    
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

%% Find outliers across all time points
% Can't do this, I toss way too many trials!
%   'TFRmth_S1t2_zS8t0_f2t40' for bhv z and nrl logz: 57, 82, 149, 64, 111, 175
%   'TFRmth_S1t2_madS8t0_f2t40' for bhv z and nrl z: 55, 110, 180, 48, 121, 174
% pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
% outlier_ix = cell(size(pow_vars));
% for t_ix = 1:numel(plt_time_vec)
%     % Normalize and combine power data per frequency band and SBJ
%     pow_data = struct;
%     for f = 1:length(pow_vars)
%         pow_data.(pow_vars{f}) = [];
%     end
%     for s = 1:length(SBJs)
%         pow_data.PFC_theta  = [pow_data.PFC_theta; fn_normalize_predictor(theta_pow{s,2}(:,t_ix),norm_nrl_pred)];
%         pow_data.PFC_betalo = [pow_data.PFC_betalo; fn_normalize_predictor(betalo_pow{s,2}(:,t_ix),norm_nrl_pred)];
%         pow_data.PFC_betahi = [pow_data.PFC_betahi; fn_normalize_predictor(betahi_pow{s,2}(:,t_ix),norm_nrl_pred)];
%         pow_data.BG_theta   = [pow_data.BG_theta; fn_normalize_predictor(theta_pow{s,1}(:,t_ix),norm_nrl_pred)];
%         pow_data.BG_betalo  = [pow_data.BG_betalo; fn_normalize_predictor(betalo_pow{s,1}(:,t_ix),norm_nrl_pred)];
%         pow_data.BG_betahi  = [pow_data.BG_betahi; fn_normalize_predictor(betahi_pow{s,1}(:,t_ix),norm_nrl_pred)];
%     end
%     
%     % Identify outliers
%     for f = 1:length(pow_vars)
%         % Identify outliers
%         out_idx = abs(pow_data.(pow_vars{f}))>outlier_thresh;
%         
%         % Track outliers per channel and frequency band
%         outlier_ix{f} = [outlier_ix{f}; find(out_idx)];
%     end
% end
% 
% % Compile list of all trials to toss per frequency band:
% outlier_ix_fband = cell(size(pow_vars));
% for f = 1:length(pow_vars)
%     outlier_ix_fband{f} = unique(outlier_ix{f});
% end

%% Find outliers tossed from time-averaged analysis
% Load group model table
table_name = [an_id norm_bhv_str norm_nrl_str];
table_all_fname = [prj_dir 'data/GRP/GRP_' table_name '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all_avg = readtable(table_all_fname);

% Identify outliers
pow_vars = {'PFC_theta','PFC_betalo','PFC_betahi','BG_theta','BG_betalo','BG_betahi'};
out_idx_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    out_idx_all.(pow_vars{f}) = abs(table_all_avg.(pow_vars{f}))>outlier_thresh;
    fprintf(2,'Bad trials in table_all for %s: %d\n',pow_vars{f},sum(out_idx_all.(pow_vars{f})));
end

%% Run model per time point
fprintf('Running LMEs: (total = %d)\n\t',length(plt_time_vec));
theta_lme = cell([2 numel(plt_time_vec)]);
beta_lme  = cell([2 numel(plt_time_vec)]);
theta_stat = cell([2 numel(plt_time_vec)]);
beta_stat  = cell([2 numel(plt_time_vec)]);
for t_ix = 1:numel(plt_time_vec)
    %% Create time-specific neural table
    PFC_theta  = [];
    PFC_betalo = [];
    PFC_betahi = [];
    BG_theta   = [];
    BG_betalo  = [];
    BG_betahi  = [];
    for s = 1:length(SBJs)
        PFC_theta  = [PFC_theta; fn_normalize_predictor(theta_pow{s,2}(:,t_ix),norm_nrl_pred)];
        PFC_betalo = [PFC_betalo; fn_normalize_predictor(betalo_pow{s,2}(:,t_ix),norm_nrl_pred)];
        PFC_betahi = [PFC_betahi; fn_normalize_predictor(betahi_pow{s,2}(:,t_ix),norm_nrl_pred)];
        BG_theta   = [BG_theta; fn_normalize_predictor(theta_pow{s,1}(:,t_ix),norm_nrl_pred)];
        BG_betalo  = [BG_betalo; fn_normalize_predictor(betalo_pow{s,1}(:,t_ix),norm_nrl_pred)];
        BG_betahi  = [BG_betahi; fn_normalize_predictor(betahi_pow{s,1}(:,t_ix),norm_nrl_pred)];
    end
    
    %% Convert into table format suitable for LME modelling
    table_all  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
        rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_diff_cur,...
        rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_diff_prv);
    
    %% Toss outliers
    good_tbl_all = struct;
    for f = 1:length(pow_vars)
        % Toss outlier trials determined from time-averaged model
        good_tbl_all.(pow_vars{f}) = table_all(~out_idx_all.(pow_vars{f}),:);
    end
    
    %% Create previous trial table (Toss NaNs)
    good_tbl_prv = good_tbl_all;
    for p = 1:length(pow_vars)
        prv_nan_idx = isnan(good_tbl_prv.(pow_vars{p}).SV_prv);
        good_tbl_prv.(pow_vars{p})(prv_nan_idx,:) = [];
        prv_fields = good_tbl_prv.(pow_vars{p}).Properties.VariableNames;
        for f = 1:length(prv_fields)
            if any(isnan(good_tbl_prv.(pow_vars{p}).(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
        end
    end
    
    %% Run LME models
    fprintf('%.02f..',plt_time_vec(t_ix));
    % PFC theta and previous reward:
    lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
    lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
    theta_stat{2,t_ix} = compare(lme0,lme1,'CheckNesting',true);%,'NSim',1000)
    theta_lme{2,t_ix} = lme1;
    
    % PFC beta low and effort:
    lme0 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ 1 + (1|sbj_n)');
    lme1 = fitlme(good_tbl_all.PFC_betalo,'PFC_betalo~ effortS_cur + (1|sbj_n)');
    beta_stat{2,t_ix} = compare(lme0,lme1,'CheckNesting',true);%,'NSim',1000)
    beta_lme{2,t_ix}  = lme1;
    
    % BG theta and previous reward:
    lme0 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
    lme1 = fitlme(good_tbl_prv.BG_theta,'BG_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
    theta_stat{1,t_ix} = compare(lme0,lme1,'CheckNesting',true);%,'NSim',1000)
    theta_lme{1,t_ix} = lme1;
    
    % BG beta low and effort:
    lme0 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ BG_roi + (1|sbj_n)');
    lme1 = fitlme(good_tbl_all.BG_betalo,'BG_betalo~ effortS_cur + BG_roi + (1|sbj_n)');
    beta_stat{1,t_ix} = compare(lme0,lme1,'CheckNesting',true);%,'NSim',1000)
    beta_lme{1,t_ix}  = lme1;
    
    if mod(t_ix,10)==0; fprintf('\n\t'); end
end
fprintf('\n');

%% Save LMEs
grp_dir = [prj_dir 'data/GRP/'];
lmm_fname = [grp_dir 'GRP_' an_id norm_bhv_str norm_nrl_str out_thresh_str '_LMM_ts.mat'];
fprintf('Saving %s\n',lmm_fname);
save(lmm_fname,'-v7.3','theta_lme','beta_lme','theta_stat','beta_stat','plt_time_vec');

