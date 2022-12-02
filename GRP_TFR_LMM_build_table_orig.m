%% Mixed modelling on the Processed data %%
% This was the version with only outliers tossed from PFC theta
%   around the time of the lab meeting on Aug 23, 2022 I adedd
%   band-specific outliers, but need to decide how to toss those trials
% this will be deleted in the next commit!
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
% an_id = 'TFRmth_S1t2_madS8t0_f2t40';%'TFRmth_S1t2_zS1t0_f2t40_log';%
an_id = 'TFRmth_D1t1_madS8t0_f2t40';
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
    elseif contains(an_id,{'TFRmth_S1t2_zS8t0_f2t40','TFRmth_S1t2_madS8t0_f2t40'})
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
    elseif contains(an_id,{'TFRmth_D1t1_zS8t0_f2t40','TFRmth_D1t1_madS8t0_f2t40'})
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
if contains(an_id,'log')
    out_thresh_str = ['out' num2str(log_outlier_thresh)];
else
    out_thresh_str = '';
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

%% Convert into table format suitable for LME modelling
table_cur  = table(trl_n_cur, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
                 rt_cur, logrt_cur, reward_cur, effort_cur, effortS_cur, decision_cur, SV_cur, absSV_cur, pAccept_cur, dec_diff_cur);
table_prv = table(trl_n_prv, sbj_n, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
                 rt_cur, logrt_cur, rt_prv, logrt_prv, reward_prv, effort_prv, effortS_prv, decision_prv, SV_prv, absSV_prv, pAccept_prv, dec_diff_prv);
prv_nan_idx = isnan(table_prv.SV_prv);
table_prv(prv_nan_idx,:) = [];
prv_fields = table_prv.Properties.VariableNames;
for f = 1:length(prv_fields)
    if any(isnan(table_prv.(prv_fields{f}))); error(['NaN is table_prv.' prv_fields{f}]); end
end

table_cur_match = table_cur;
table_cur_match(prv_nan_idx,:) = [];
table_cur_match.sbj_n      = [];
table_cur_match.rt_cur     = [];
table_cur_match.logrt_cur  = [];
table_cur_match.PFC_roi    = [];
table_cur_match.BG_roi     = [];
table_cur_match.PFC_theta  = [];
table_cur_match.PFC_betalo = [];
table_cur_match.PFC_betahi = [];
table_cur_match.BG_theta   = [];
table_cur_match.BG_betalo  = [];
table_cur_match.BG_betahi  = [];
table_all = [table_cur_match, table_prv];

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
        toss_ix = find(abs(table_all.PFC_theta)>log_outlier_thresh);
        table_all(toss_ix,:)  = [];
    else
        fprintf('No bad trials for PFC theta with threshold %f\n',log_outlier_thresh);
    end
end

% DatTable3 = [table_cur, DatTable2tmp];
% table_all = [table_cur, table_prv_only];
% % These are different by trials 1 and 2 for each patient, unclear why (they
% % look identical...)

%% Write tables for R
% current trial table
table_cur_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_cur.csv'];
fprintf('\tSaving %s...\n',table_cur_fname);
writetable(table_cur,table_cur_fname);

% previous trial table
table_prv_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_prv.csv'];
fprintf('\tSaving %s...\n',table_prv_fname);
writetable(table_prv,table_prv_fname);

% combined table
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id out_thresh_str norm_bhv_str norm_nrl_str '_full_table_all.csv'];
fprintf('\tSaving %s...\n',table_all_fname);
writetable(table_all,table_all_fname);

