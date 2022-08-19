%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%%
SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};
man_trl_rej_ix = {[], [71 72], [], [27 28 79 80 86 87 97 98 102 103 128 139 140 148 149 150]};
% sbj_colors = distinguishable_colors(length(SBJs));

% Analysis parameters:
norm_bhv_reg = 0;
an_id = 'simon_S';%'TFRw_S25t2_dbS25t05_fl2t40_c7';%'TFRw_D1t1_dbS25t05_fl2t40_c7';%
toss_same_trials = 1;
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
    plt_lim = [-0.2 2];
    betahi_cf = ones([2 numel(SBJs)])*-1;
    if contains(an_id,'simon')
        % Simon params: [10,17,13,12]
        theta_cf = [-1 -1 -1 -1; -1 -1 -1 -1]; % BG (row1) then PFC (row2)
        betalo_cf = [-1 -1 -1 -1; 10 17 13 12]; % PFC03, PFC04, PFC05, PFC01
    elseif strcmp(an_id,'TFRw_S25t2_dbS25t05_fl2t40_c7')
        theta_cf = [2.5 3.5 3.5 6; 5 4 3.5 3]; % BG (row1) then PFC (row2)
        betalo_cf  = [-1 -1 -1 -1; -1 -1 -1 -1];
    else
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
        theta_cf = [2.5 3 3.5 5; 5 3.5 3.5 3]; % BG (row1) then PFC (row2)
        betalo_cf  = [17 22 14 14; 13 17 11 15];
        %     betahi_cf  = % no SBJ-specific peaks;
        %     betalo_cf = ones([2 numel(SBJs)])*-1;
    end
elseif contains(an_id,'_D')
    an_lim = [-0.25 0.25];
    error('not ready, simon results have no beta desynch');
    if strcmp(an_id,'TFRw_D1t1_dbS25t05_fl2t40_c7')
        theta_cf = [3.5 3 4.5 6; 3 3 4.5 3]; % BG (row1) then PFC (row2)
        betalo_cf  = [13 20 17 18; -1 -1 -1 -1];
        betahi_cf = ones([2 numel(SBJs)])*-1;
    else
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
        theta_cf = [3.5 3 3.5 7; 4.5 3 3.5 3]; % BG (row1) then PFC (row2)
        %     betalo_cf  = [17 22 14 14; 13 17 11 15];
        %     betahi_cf  = % no SBJ-specific peaks;
        betalo_cf = [20 17 11 24; -1 -1 -1 -1];
        betahi_cf = ones([2 numel(SBJs)])*-1;
    end
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
if contains(an_id,'simon')
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
    if contains(an_id,'simon')
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
    
    %% Initialize data
    full_time_vec = tfr.time;
    full_freq_vec = tfr.freq;
    
    % Check channel index
    if ~strcmp(tfr.label{1},'LFP'); error('BG LFP is not first channel!'); end
    if ~any(strcmp(tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    %% Compute single trial power
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
            fprintf(2,'\tWARNING: %s theta center freq is off: aim = %.02f, actual = %.02f...\n',...
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

%% Create data table for behavior
% Initialise variables (A = current trial, As = previous trial)
sbj_n_A = [];
PFC_roi = [];
BG_roi  = [];
trl_n_A = [];

% StakeA=[];
reward_A   = [];
effort_A   = [];
EFF_A      = [];
decision_A = [];
SV_A       = [];
SV_md_A    = [];

reward_As   = [];
effort_As   = [];
effortO_As  = [];
decision_As = [];
SV_As       = [];

for s = 1:4
    % Concatenate SBJ, beta, theta values
    trl_n = size(bhvs{s}.stake,1);
    trl_n_A = [trl_n_A; [1:trl_n]'];
    sbj_n_A = [sbj_n_A; num2str(ones(trl_n,1).*s)];
    if strcmp(sbj_pfc_roi{s},'OFC'); pfc_roi_ix = 1; else; pfc_roi_ix = 2; end
    PFC_roi = [PFC_roi; num2str(ones(trl_n,1).*pfc_roi_ix)];
    if strcmp(sbj_bg_roi{s},'STN'); bg_roi_ix = 1; else; bg_roi_ix = 2; end
    BG_roi = [BG_roi; num2str(ones(trl_n,1).*bg_roi_ix)];
    
    % Extract the behavioral variables %
    reward_A   = [reward_A; bhvs{s}.stake];
    effort_A   = [effort_A; bhvs{s}.effort];
    EFF_A      = [EFF_A; bhvs{s}.EFFs];
    decision_A = [decision_A; bhvs{s}.decision];
    SV_A       = [SV_A; bhvs{s}.SV];
    SV_md_A    = [SV_md_A; abs(bhvs{s}.SV-median(bhvs{s}.SV))]; % maybe conflict? distance from indecision point? difficulty?
    
    % Reward
    rew_prv = smooth(bhvs{s}.stake,2); % why isn't this doing anything?
    rew_prv(end-1:end)=[]; rew_prv=[nan(2,1); rew_prv]; % why shift by 2 instead of 1?
    reward_As = [reward_As; rew_prv];%bhvs{s}.stake_prv];
    
    % Effort
    effort_prv = smooth(bhvs{s}.effort,1);
    effort_prv(end)=[]; effort_prv = [nan(1,1); effort_prv];
    effort_As = [effort_As; effort_prv];%bhvs{s}.effort_prv];
    
    % Obj. Effort
    %     effortO_prv = smooth(bhvs{s}.EFFs,1);
    effortO_prv = bhvs{s}.EFFs;
    effortO_prv(end)=[]; effortO_prv = [nan(1,1); effortO_prv];
    effortO_As = [effortO_As; effortO_prv];
    
    % decision_A
    dec_prv = smooth(bhvs{s}.decision,1);
    dec_prv(end) = []; dec_prv = [nan(1,1); dec_prv];
    decision_As = [decision_As; dec_prv];%bhvs{s}.decision_prv];
    
    % Subj Val.
    %     SV_prv = smooth(bhvs{s}.SV,1);
    SV_prv = bhvs{s}.SV;
    SV_prv(end) = []; SV_prv = [nan(1,1); SV_prv];
    SV_As = [SV_As; SV_prv];
end

%% Run model per time point
fprintf('Running LMEs:\n\t');
theta_lme = cell([2 numel(plt_time_vec)]);
beta_lme  = cell([2 numel(plt_time_vec)]);
for t_ix = 1:numel(plt_time_vec)
    %% Convert into table format suitable for LME modelling
    % Initialise variables (A = current trial, As = previous trial)
    PFC_theta = [];
    PFC_betalo  = [];
    PFC_betahi  = [];
    BG_theta  = [];
    BG_betalo   = [];
    BG_betahi   = [];
    
    % Create time-specific neural table
    for s = 1:length(SBJs)
        PFC_theta  = [PFC_theta; theta_pow{s,2}(:,t_ix)];
        PFC_betalo = [PFC_betalo; betalo_pow{s,2}(:,t_ix)];
        PFC_betahi = [PFC_betahi; betahi_pow{s,2}(:,t_ix)];
        BG_theta   = [BG_theta; theta_pow{s,1}(:,t_ix)];
        BG_betalo  = [BG_betalo; betalo_pow{s,1}(:,t_ix)];
        BG_betahi  = [BG_betahi; betahi_pow{s,1}(:,t_ix)];
    end
    
    table_A  = table(trl_n_A, sbj_n_A, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
        reward_A, effort_A, decision_A, SV_A, EFF_A);
    table_As = table(trl_n_A, sbj_n_A, PFC_roi, BG_roi, PFC_theta, PFC_betalo, PFC_betahi, BG_theta, BG_betalo, BG_betahi,...
        reward_As, effort_As, decision_As, SV_As, effortO_As);
    
    table_As_only = table(reward_As, effort_As, decision_As, SV_As, effortO_As);
    
    DatTable2tmp = table_As;
    DatTable2tmp.sbj_n_A    = [];
    DatTable2tmp.trl_n_A    = [];
    DatTable2tmp.PFC_roi    = [];
    DatTable2tmp.BG_roi     = [];
    DatTable2tmp.PFC_theta  = [];
    DatTable2tmp.PFC_betalo = [];
    DatTable2tmp.PFC_betahi = [];
    DatTable2tmp.BG_theta   = [];
    DatTable2tmp.BG_betalo  = [];
    DatTable2tmp.BG_betahi  = [];
    
    DatTable3 = [table_A, DatTable2tmp];
    table_all = [table_A, table_As_only];
    % These are different by trials 1 and 2 for each patient, unclear why (they
    % look identical...)
    
    %% Run LME models
    fprintf('%.02f..',plt_time_vec(t_ix));
    theta_lme{1,t_ix} = fitlme(table_As,'BG_theta~ SV_As*PFC_roi + (1|sbj_n_A)');
    beta_lme{1,t_ix}  = fitlme(table_A,'BG_betalo~ SV_A*PFC_roi + (1|sbj_n_A)');
    theta_lme{2,t_ix} = fitlme(table_As,'PFC_theta~ SV_As*PFC_roi + (1|sbj_n_A)');
    beta_lme{2,t_ix}  = fitlme(table_A,'PFC_betalo~ SV_A*PFC_roi + (1|sbj_n_A)');
    if mod(t_ix,10)==0; fprintf('\n\t'); end
end
fprintf('\n');

%% Save LMEs
grp_dir = [prj_dir 'data/GRP/'];
lme_fname = [grp_dir 'GRP_' an_id '_LME_ts.mat'];
fprintf('Saving %s\n',lme_fname);
save(lme_fname,'-v7.3','theta_lme','beta_lme','plt_time_vec');

%% New LME Modeling
% % Comparisons: effort and reward separately
% 
% % ================ PFC ================
% % Current trial Subjective Value
% lme = fitlme(table_A,'PFC_theta~ SV_A*PFC_roi + (1|sbj_n_A)')   % PFC_roi p = 0.02; D -0.25:0.25 SV is trending p=0.056
% lme = fitlme(table_A,'PFC_betalo~ SV_A*PFC_roi + (1|sbj_n_A)')
% lme = fitlme(table_A,'PFC_betahi~ SV_A*PFC_roi + (1|sbj_n_A)')  % D -0.25:0.25 SV p = 0.044
% % lme = fitlme(table_A,'PFC_betaHi~ SV_A*PFC_roi + (1|sbj_n_A)')
% 
% % Previous trial Subjective Value
% lme = fitlme(table_As,'PFC_theta~ SV_As*PFC_roi + (1|sbj_n_A)') % PFC_roi p = 0.0006, PFC*SV p = 0.016
% lme = fitlme(table_As,'PFC_betalo~ SV_As*PFC_roi + (1|sbj_n_A)')
% lme = fitlme(table_As,'PFC_betahi~ SV_As*PFC_roi + (1|sbj_n_A)')
% 
% % ================ BG ================
% % Current trial Subjective Value
% lme = fitlme(table_A,'BG_theta~ SV_A*BG_roi + (1|sbj_n_A)')     % BG_roi p = 0.03
% lme = fitlme(table_A,'BG_betalo~ SV_A*BG_roi + (1|sbj_n_A)')    % BG_roi p = 0.0004; D -0.25:0.25 SV p - 0.03
% lme = fitlme(table_A,'BG_betahi~ SV_A*BG_roi + (1|sbj_n_A)')
% 
% % Previous trial Subjective Value
% lme = fitlme(table_As,'BG_theta~ SV_As*BG_roi + (1|sbj_n_A)')   % BG_roi p = 0.04
% lme = fitlme(table_As,'BG_betalo~ SV_As*BG_roi + (1|sbj_n_A)')  % BG_roi p = 0.0009; D -0.25:0.25 SV_As p = 0.01
% lme = fitlme(table_As,'BG_betahi~ SV_As*BG_roi + (1|sbj_n_A)')
% 
% %% LME Modelling %%
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
