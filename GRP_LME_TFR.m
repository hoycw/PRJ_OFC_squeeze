%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};
man_trl_rej_ix = {[], [71 72], [], [27 28 79 80 86 87 97 98 102 103 128 139 140 148 149 150]};
% sbj_colors = distinguishable_colors(length(SBJs));

% Analysis parameters:
norm_bhv_pred = 'zscore';%'none';%
norm_nrl_pred = 'zscore';%'none';%
log_outlier_thresh = 7;
an_id = 'TFRmth_S1t2_zS8t0_f2t40_log';%'TFRmth_S1t2_zS1t0_f2t40_log';%
% an_id = 'TFRmth_D1t1_zS8t0_f2t40';
use_simon_tfr = 0;
toss_same_trials = 1;
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
    betahi_cf = ones([2 numel(SBJs)])*-1;
    if use_simon_tfr
        % Simon params: [10,17,13,12]
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf = [-1 -1 -1 -1; 10 17 13 12]; % PFC03, PFC04, PFC05, PFC01
    elseif strcmp(an_id,'TFRmth_S1t2_zS1t0_f2t40')  % equivalent to Simon
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; 10 17 13 12];
    elseif contains(an_id,'TFRmth_S1t2_zS8t0_f2t40')
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
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
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

%% Convert into table format suitable for LME modelling
% Initialise variables (cur = current trial, prv = previous trial)
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

table_cur  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
                 rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_diff_cur);
table_prv = table(trl_n_prv, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
                 rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_diff_prv);
table_prv(isnan(table_prv.SV_prv),:) = [];
prv_fields = table_prv.Properties.VariableNames;
for f = 1:length(prv_fields)
    if any(isnan(table_prv.(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
end
% table_prv_only = table(reward_prv, effort_prv, decision_prv, SV_prv, effortS_prv, rt_prv, logrt_prv);

% Toss theta outlier
if contains(an_id,'log')
    if any(abs(table_prv.PFC_theta)>log_outlier_thresh)
        toss_ix = find(abs(table_prv.PFC_theta)>log_outlier_thresh);
        fprintf(2,'\tTossing trials for PFC theta:');
        disp(table_prv.PFC_theta(toss_ix));
        fprintf('\n');
        table_prv(toss_ix,:) = [];
        toss_ix = find(abs(table_cur.PFC_theta)>log_outlier_thresh);
        table_cur(toss_ix,:)  = [];
    else
        fprintf('No bad trials for PFC theta with threshold %f\n',log_outlier_thresh);
    end
end

% DatTable2tmp = table_prv;
% DatTable2tmp.sbj_n    = [];
% DatTable2tmp.trl_n_cur    = [];
% DatTable2tmp.PFC_roi    = [];
% DatTable2tmp.BG_roi     = [];
% DatTable2tmp.PFC_theta  = [];
% DatTable2tmp.PFC_betalo = [];
% DatTable2tmp.PFC_betahi = [];
% DatTable2tmp.BG_theta   = [];
% DatTable2tmp.BG_betalo  = [];
% DatTable2tmp.BG_betahi  = [];
% 
% DatTable3 = [table_cur, DatTable2tmp];
% table_all = [table_cur, table_prv_only];
% % These are different by trials 1 and 2 for each patient, unclear why (they
% % look identical...)

%% Write tables for R
% current trial table
table_cur_fname = [prj_dir 'data/GRP/GRP_' an_id norm_bhv_str norm_nrl_str '_full_table_cur.csv'];
fprintf('\tSaving %s...\n',table_cur_fname);
writetable(table_cur,table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' an_id norm_bhv_str norm_nrl_str '_full_table_prv.csv'];
fprintf('\tSaving %s...\n',table_prv_fname);
writetable(table_prv,table_prv_fname);

%% LME Model Selection
% ================ PFC ================
% Current trial Subjective Value
% lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1|sbj_n)');
% lme2 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1+SV_cur|sbj_n)');
% lme3 = fitlme(table_cur,'PFC_betalo~ SV_cur + PFC_roi + (1|sbj_n)');
% lme4 = fitlme(table_cur,'PFC_betalo~ SV_cur*PFC_roi + (1|sbj_n)');
% pfc_beta_sbj_slope = compare(lme1,lme2,'NSim',1000);
% pfc_beta_roi_fix   = compare(lme1,lme3,'NSim',1000);
% pfc_beta_roi_int   = compare(lme3,lme4,'NSim',1000);
% pfc_beta_sv_fix    = compare(lme0,lme1,'NSim',1000)

% Previous trial Subjective Value
% lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');
% lme1 = fitlme(table_prv,'PFC_theta~ SV_prv + (1|sbj_n)');
% lme2 = fitlme(table_prv,'PFC_theta~ SV_prv + (1+SV_prv|sbj_n)');
% lme3 = fitlme(table_prv,'PFC_theta~ SV_prv + PFC_roi + (1|sbj_n)');
% lme4 = fitlme(table_prv,'PFC_theta~ SV_prv*PFC_roi + (1|sbj_n)');
% pfc_theta_sbj_slope = compare(lme1,lme2,'NSim',1000);
% pfc_theta_roi_fix   = compare(lme1,lme3,'NSim',1000);
% pfc_theta_roi_int   = compare(lme3,lme4,'NSim',1000);
% pfc_theta_svp_fix   = compare(lme0,lme1,'NSim',1000)

% ================ BG ================
% Current trial Subjective Value
% lme0 = fitlme(table_cur,'BG_betalo~ BG_roi + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(table_cur,'BG_betalo~ SV_cur + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(table_cur,'BG_betalo~ SV_cur + (1+SV_cur|sbj_n)');%,'StartMethod','random');
% lme3 = fitlme(table_cur,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)');%,'StartMethod','random');
% lme4 = fitlme(table_cur,'BG_betalo~ SV_cur*BG_roi + (1|sbj_n)');%,'StartMethod','random');
% bg_beta_sbj_slope = compare(lme1,lme2,'NSim',1000);
% bg_beta_roi_fix   = compare(lme1,lme3,'NSim',1000);
% bg_beta_roi_int   = compare(lme3,lme4,'NSim',1000);
% bg_beta_sv_fix    = compare(lme0,lme3,'NSim',1000)
% roi_fix_pvals = nan([100 1]);
% for m = 1:100
%     lme1 = fitlme(table_cur,'BG_betalo~ SV_cur + (1|sbj_n)','StartMethod','random');
%     lme3 = fitlme(table_cur,'BG_betalo~ SV_cur + BG_roi + (1|sbj_n)','StartMethod','random');
%     roi_fix   = compare(lme1,lme3,'NSim',1000);
%     roi_fix_pvals(m) = roi_fix.pValue(2);
% end

% Previous trial Subjective Value
% lme0 = fitlme(table_prv,'BG_theta~ BG_roi + (1|sbj_n)');%,'StartMethod','random');
% lme1 = fitlme(table_prv,'BG_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
% lme2 = fitlme(table_prv,'BG_theta~ SV_prv + (1+SV_prv|sbj_n)');%,'StartMethod','random');
% lme3 = fitlme(table_prv,'BG_theta~ SV_prv + BG_roi + (1|sbj_n)');%,'StartMethod','random');
% lme4 = fitlme(table_prv,'BG_theta~ SV_prv*BG_roi + (1|sbj_n)');%,'StartMethod','random');
% bg_theta_sbj_slope = compare(lme1,lme2,'NSim',1000);
% bg_theta_roi_fix   = compare(lme1,lme3,'NSim',1000);
% bg_theta_roi_int   = compare(lme3,lme4,'NSim',1000);
% bg_theta_svp_fix   = compare(lme0,lme3,'NSim',1000)

%% Test best models:
% PFC beta low and subejctive value:
lme0 = fitlme(table_cur,'PFC_betalo~ 1 + (1|sbj_n)');
lme1 = fitlme(table_cur,'PFC_betalo~ SV_cur + (1|sbj_n)');
pfc_betalo_sv = compare(lme0,lme1)%,'NSim',1000)

% PFC theta and previous subjective value:
lme0 = fitlme(table_prv,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
lme1 = fitlme(table_prv,'PFC_theta~ SV_prv + (1|sbj_n)');%,'StartMethod','random');
pfc_theta_sv_prv = compare(lme0,lme1)%,'NSim',1000)

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


%% LME Modelling %%
% % OFC:
% % Current trial
% lme = fitlme(table_A,'PFC_betalo~ SV_A + (1|sbj_n_A)')     %SV p = 0.01
% lme = fitlme(table_A,'PFC_theta~ SV_A + (1|sbj_n_A)')    %SV p = 0.4
% 
% lme = fitlme(table_A,'PFC_betalo~ SV_A*PFC_roi + (1|sbj_n_A)')     %SV p = 0.0007; pfc_roi p = 0.0048; * p = 0.1
% % lme = fitlme(table_A,'PFC_beta~ SV_A + (SV_A|pfc_roi) + (1|sbj_n_A)')     %SV p = 0.048
% lme = fitlme(table_A,'PFC_theta~ SV_A*PFC_roi + (1|sbj_n_A)')    %none
% 
% return;
% % Previous trial
% lme = fitlme(table_As,'PFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04
% lme = fitlme(table_As,'PFC_beta~ SV_As + (1|sbj_n_A)')   %
% lme = fitlme(table_As,'BG_theta~ SV_As + (1|sbj_n_A)')
% lme = fitlme(table_As,'BG_beta~ SV_As + (1|sbj_n_A)')
% 
% lme = fitlme(table_As,'PFC_theta~ SV_As*pfc_roi + (1|sbj_n_A)')  %SV p = 0.052
% lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|pfc_roi) + (1|sbj_n_A)')  %SV p = 0.038
% 
% lme = fitlme(table_A,'BG_beta~ SV_A + (1|sbj_n_A)')
% lme = fitlme(table_A,'BG_theta~ SV_A + (1|sbj_n_A)')
% 
% 
% lme = fitlme(table_As,'PFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04
% 
% 
% % predict decision from neural data
% glme = fitglme(DatTable3, 'decision_A~ PFC_beta + PFC_theta + (1|sbj_n_A)','Link','log')    % beta p = 0.017
% glme = fitglme(DatTable3, 'decision_A~ BG_beta + BG_theta + (1|sbj_n_A)','Link','log')
% 
% 
% return;
% 
% lme = fitlme(DatTable,'PFC_theta~ SV_A + (SV_A|sbj_n_A)')
% lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A)')
% lme = fitlme(DatTable,'PFC_beta~ SV_A*decision_A + (SV_A|sbj_n_A) + (decision_A|sbj_n_A)')
% lme = fitlme(table_As,'PFC_theta~ SV_As*decision_As + (SV_As|sbj_n_A) + (decision_As|sbj_n_A)')
% lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|sbj_n_A)')
% lme = fitlme(table_As,'PFC_theta~ SV_As*decision_As + (1|sbj_n_A) + (1|sbj_n_A)')
% lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|sbj_n_A)')
% lme = fitlme(table_As,'PFC_theta~ decision_As + (decision_As|sbj_n_A)')
% 
% 
% % lme = fitlme(table_As,'PFC_theta~ decision_As + (decision_As|sbj_n_A)')
% % 
% 
% 
% return;
% 
% % Objective value%
% % ObjVal=RewardA-EffortA
% 
% % lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A)')
% % lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A) + (SV_A-1|sbj_n_A)')
% % figure();
% % plotResiduals(lme,'fitted')
% % find(residuals(lme) > 1.5)
% 
% lme2 = fitlme(DatTable,'PFC_beta ~ SV_A ')
% compare(lme2,lme)
% 
% [~,~,stats]=covarianceParameters(lme)
% 
% %% Does theta = conflict?
% 
% % Rescale reward and effort by max and minimum %
% maxR=max(RewardA); minR=min(RewardA);
% maxE=max(EffortA.^2); minE=min(EffortA.^2);
% 
% maxRr=repmat(maxR,size(RewardA,1), size(RewardA,2));
% RewardArs=(RewardA./maxRr).*100;
% 
% maxEr=repmat(maxE,size(EffortA,1), size(EffortA,2));
% EffortArs=((EffortA.^2)./maxEr).*100;
% 
% % Conflict 
% % Abs here is critical - absolute conflict - distinguishes from subjective
% % value. 
% Conflict=abs(RewardArs-EffortArs);
% % Add new column to datTable
% DatTable_c=DatTable;
% DatTable_c.Conflict = Conflict;
% lmec = fitlme(DatTable_c,'PFC_theta~ Conflict + (Conflict|sbj_n_A)')
% % Conflict here is going to be related to SV - high reward 
% 
% OV=RewardArs-EffortArs;
% DatTable_co=DatTable_c;
% DatTable_co.OV = OV;
% lmeo = fitlme(DatTable_co,'PFC_beta~ OV + (OV|sbj_n_A)'); %OV p = 0.02
% 
% clc; close all;
% DatTable_cod=DatTable_co;
% DatTable_cod.DiffA = DiffA;
% lmdiff = fitlme(DatTable_cod,'PFC_theta~ DiffA + (DiffA|ptnumA)')
