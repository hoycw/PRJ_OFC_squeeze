%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%%
SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};
sbj_colors = distinguishable_colors(length(SBJs));

% Analysis parameters:
norm_bhv_reg = 0;
an_id = 'TFRw_S25t2_zbtS25t05_fl2t40_c7';%'TFRw_D1t1_zbtS25t05_fl2t40_c7';%
if contains(an_id,'_S')
    an_lim = [0.5 1.5];
elseif contains(an_id,'_D')
    an_lim = [-0.5 0];
end

thetapk_lfp = [3.5 8 3.5 5];
thetapk_pfc = [5 3.5 3.5 4.5];
betapk_lfp  = [17 22 14 15];
betapk_pfc  = [14 17 14 14];
% theta_lim  = [4 7];
% beta_lim   = [13 30];    
% simon_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
thetapk_bw = 4;
betapk_bw = 4;

%% Analysis Set Up
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

thetapk_lim = nan(length(SBJs),2,2);
betapk_lim  = nan(length(SBJs),2,2);
for s = 1:length(SBJs)
    thetapk_lim(s,1,:) = [thetapk_lfp(s)-thetapk_bw/2 thetapk_lfp(s)+thetapk_bw/2];
    thetapk_lim(s,2,:) = [thetapk_pfc(s)-thetapk_bw/2 thetapk_pfc(s)+thetapk_bw/2];
    betapk_lim(s,1,:)  = [betapk_lfp(s)-betapk_bw/2 betapk_lfp(s)+betapk_bw/2];
    betapk_lim(s,2,:)  = [betapk_pfc(s)-betapk_bw/2 betapk_pfc(s)+betapk_bw/2];
end

% Plotting parameters
sem_alpha  = 0.5;
time_lim = [0.5 1.5];
% stim_time_lim = [0 2];
% plot_freq_lim = [2 30];

save_fig = 0;

%% Compute Theta and Beta power
thetapk_pow = cell([numel(SBJs) 2]);
betapk_pow  = cell([numel(SBJs) 2]);
bhvs        = cell(size(SBJs));
for s = 1:length(SBJs)
    %% Load TFR
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_' an_id '.mat'];
    fprintf('Loading %s\n',proc_fname);
    tmp = load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    bhvs{s} = sbj_data.bhv;
    
    % Initialize data
    time_vec = tmp.tfr.time;
    freq_vec = tmp.tfr.freq;
    
    % Check channel index
    if ~strcmp(tmp.tfr.label{1},'LFP'); error('BG LFP is not first channel!'); end
    if ~any(strcmp(tmp.tfr.label{2},{'FPC','OFC'})); error('PFC is not second channel!'); end
    
    % Compute single trial power
    for ch_ix = 1:2
        % Theta
        cfg = [];
        cfg.channel = tmp.tfr.label(ch_ix);
        cfg.avgoverfreq = 'yes';
        cfg.avgovertime = 'yes';
        cfg.avgoverrpt  = 'no';
        cfg.latency     = an_lim;
        cfg.frequency   = squeeze(thetapk_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tmp.tfr);
        thetapk_pow{s,ch_ix} = pow.powspctrm;
        
        cfg.frequency = squeeze(betapk_lim(s,ch_ix,:))';
        pow = ft_selectdata(cfg, tmp.tfr);
        betapk_pow{s,ch_ix} = pow.powspctrm;
    end
end

%% Convert into table format suitable for LME modelling
% Initialise variables (A = current trial, As = previous trial)
sbj_n_A = [];
PFC_roi = [];
BG_roi  = [];
trl_n_A = [];
PFC_theta = [];
PFC_beta  = [];
BG_theta  = [];
BG_beta   = [];

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

for s = 1:length(SBJs)
    % Concatenate SBJ, beta, theta values
    trl_n = size(bhvs{s}.trl,1);
    trl_n_A = [trl_n_A; [1:trl_n]'];
    sbj_n_A = [sbj_n_A; num2str(ones(trl_n,1).*s)];
    if strcmp(sbj_pfc_roi{s},'OFC'); pfc_roi_ix = 1; else; pfc_roi_ix = 2; end
    PFC_roi = [PFC_roi; num2str(ones(trl_n,1).*pfc_roi_ix)];
    if strcmp(sbj_bg_roi{s},'STN'); bg_roi_ix = 1; else; bg_roi_ix = 2; end
    BG_roi = [BG_roi; num2str(ones(trl_n,1).*bg_roi_ix)];
    
    PFC_theta = [PFC_theta; thetapk_pow{s,2}];
    PFC_beta  = [PFC_beta; betapk_pow{s,2}];
    BG_theta  = [BG_theta; thetapk_pow{s,1}];
    BG_beta   = [BG_beta; betapk_pow{s,1}];
    
    % Extract the behavioral variables %
    reward_A   = [reward_A; bhvs{s}.stake];
    effort_A   = [effort_A; bhvs{s}.effort];
    EFF_A      = [EFF_A; bhvs{s}.EFFs];
    decision_A = [decision_A; bhvs{s}.decision];
    SV_A       = [SV_A; bhvs{s}.SV];
    SV_md_A    = [SV_md_A; abs(bhvs{s}.SV-median(bhvs{s}.SV))];
    
    % Reward
%     rew_prv = smooth(bhvs{s}.stake,2); % why isn't this doing anything?
%     rew_prv(end-1:end)=[]; rew_prv=[nan(2,1); rew_prv]; % why shift by 2 instead of 1?
    reward_As = [reward_As; bhvs{s}.stake_prv];
    
    % Effort
%     effort_prv = smooth(bhvs{s}.effort,1);
%     effort_prv(end)=[]; effort_prv = [nan(1,1); effort_prv];
    effort_As = [effort_As; bhvs{s}.effort_prv];
    
    % Obj. Effort
%     effortO_prv = smooth(bhvs{s}.EFFs,1);
    effortO_prv = bhvs{s}.EFFs;
    effortO_prv(end)=[]; effortO_prv = [nan(1,1); effortO_prv];
    effortO_As = [effortO_As; effortO_prv];
    
    % decision_A
%     dec_prv = smooth(bhvs{s}.decision,1);
%     dec_prv(end) = []; dec_prv = [nan(1,1); dec_prv];
    decision_As = [decision_As; bhvs{s}.decision_prv];
    
    % Subj Val.
%     SV_prv = smooth(bhvs{s}.SV,1);
    SV_prv = bhvs{s}.SV;
    SV_prv(end) = []; SV_prv = [nan(1,1); SV_prv];
    SV_As = [SV_As; SV_prv];
end

table_A  = table(trl_n_A, sbj_n_A, PFC_roi, BG_roi, PFC_theta, PFC_beta, BG_theta, BG_beta,...
                 reward_A, effort_A, decision_A, SV_A, EFF_A);
table_As = table(trl_n_A, sbj_n_A, PFC_roi, BG_roi, PFC_theta, PFC_beta, BG_theta, BG_beta,...
                 reward_As, effort_As, decision_As, SV_As, effortO_As);

table_As_only = table(reward_As, effort_As, decision_As, SV_As, effortO_As);

DatTable2tmp = table_As;
DatTable2tmp.sbj_n_A   = [];
DatTable2tmp.trl_n_A   = [];
DatTable2tmp.PFC_roi   = [];
DatTable2tmp.BG_roi    = [];
DatTable2tmp.PFC_theta = [];
DatTable2tmp.PFC_beta  = [];
DatTable2tmp.BG_theta  = [];
DatTable2tmp.BG_beta   = [];

DatTable3 = [table_A, DatTable2tmp];
table_all = [table_A, table_As_only];
% These are different by trials 1 and 2 for each patient, unclear why (they
% look identical...)

%% New LME Modeling
% Comparisons: effort and reward separately
% Add high beta

% ================ PFC ================
% Current trial Subjective Value
lme = fitlme(table_A,'PFC_beta~ SV_A*PFC_roi + (1|sbj_n_A)');
lme = fitlme(table_A,'PFC_theta~ SV_A*PFC_roi + (1|sbj_n_A)')
% lme = fitlme(table_A,'PFC_betaHi~ SV_A*PFC_roi + (1|sbj_n_A)')

% Previous trial Subjective Value
lme = fitlme(table_As,'PFC_theta~ SV_As*PFC_roi + (1|sbj_n_A)')
lme = fitlme(table_As,'PFC_beta~ SV_As*PFC_roi + (1|sbj_n_A)')

% ================ BG ================
% Current trial Subjective Value
lme = fitlme(table_A,'BG_beta~ SV_A*BG_roi + (1|sbj_n_A)')
lme = fitlme(table_A,'BG_theta~ SV_A*BG_roi + (1|sbj_n_A)')

% Previous trial Subjective Value
lme = fitlme(table_As,'BG_theta~ SV_As*BG_roi + (1|sbj_n_A)')
lme = fitlme(table_As,'BG_beta~ SV_As*BG_roi + (1|sbj_n_A)')

%% LME Modelling %%
% OFC:
% Current trial
lme = fitlme(table_A,'PFC_beta~ SV_A + (1|sbj_n_A)')     %SV p = 0.01
lme = fitlme(table_A,'PFC_theta~ SV_A + (1|sbj_n_A)')    %SV p = 0.4

lme = fitlme(table_A,'PFC_beta~ SV_A*pfc_roi + (1|sbj_n_A)')     %SV p = 0.0007; pfc_roi p = 0.0048; * p = 0.1
% lme = fitlme(table_A,'PFC_beta~ SV_A + (SV_A|pfc_roi) + (1|sbj_n_A)')     %SV p = 0.048
lme = fitlme(table_A,'PFC_theta~ SV_A*pfc_roi + (1|sbj_n_A)')    %none

return;
% Previous trial
lme = fitlme(table_As,'PFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04
lme = fitlme(table_As,'PFC_beta~ SV_As + (1|sbj_n_A)')   %
lme = fitlme(table_As,'BG_theta~ SV_As + (1|sbj_n_A)')
lme = fitlme(table_As,'BG_beta~ SV_As + (1|sbj_n_A)')

lme = fitlme(table_As,'PFC_theta~ SV_As*pfc_roi + (1|sbj_n_A)')  %SV p = 0.052
lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|pfc_roi) + (1|sbj_n_A)')  %SV p = 0.038

lme = fitlme(table_A,'BG_beta~ SV_A + (1|sbj_n_A)')
lme = fitlme(table_A,'BG_theta~ SV_A + (1|sbj_n_A)')


lme = fitlme(table_As,'PFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04


% predict decision from neural data
glme = fitglme(DatTable3, 'decision_A~ PFC_beta + PFC_theta + (1|sbj_n_A)','Link','log')    % beta p = 0.017
glme = fitglme(DatTable3, 'decision_A~ BG_beta + BG_theta + (1|sbj_n_A)','Link','log')


return;

lme = fitlme(DatTable,'PFC_theta~ SV_A + (SV_A|sbj_n_A)')
lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A)')
lme = fitlme(DatTable,'PFC_beta~ SV_A*decision_A + (SV_A|sbj_n_A) + (decision_A|sbj_n_A)')
lme = fitlme(table_As,'PFC_theta~ SV_As*decision_As + (SV_As|sbj_n_A) + (decision_As|sbj_n_A)')
lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|sbj_n_A)')
lme = fitlme(table_As,'PFC_theta~ SV_As*decision_As + (1|sbj_n_A) + (1|sbj_n_A)')
lme = fitlme(table_As,'PFC_theta~ SV_As + (SV_As|sbj_n_A)')
lme = fitlme(table_As,'PFC_theta~ decision_As + (decision_As|sbj_n_A)')


% lme = fitlme(table_As,'PFC_theta~ decision_As + (decision_As|sbj_n_A)')
% 


return;

% Objective value%
% ObjVal=RewardA-EffortA

% lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A)')
% lme = fitlme(DatTable,'PFC_beta~ SV_A + (1|sbj_n_A) + (SV_A-1|sbj_n_A)')
% figure();
% plotResiduals(lme,'fitted')
% find(residuals(lme) > 1.5)

lme2 = fitlme(DatTable,'PFC_beta ~ SV_A ')
compare(lme2,lme)

[~,~,stats]=covarianceParameters(lme)

%% Does theta = conflict?

% Rescale reward and effort by max and minimum %
maxR=max(RewardA); minR=min(RewardA);
maxE=max(EffortA.^2); minE=min(EffortA.^2);

maxRr=repmat(maxR,size(RewardA,1), size(RewardA,2));
RewardArs=(RewardA./maxRr).*100;

maxEr=repmat(maxE,size(EffortA,1), size(EffortA,2));
EffortArs=((EffortA.^2)./maxEr).*100;

% Conflict 
% Abs here is critical - absolute conflict - distinguishes from subjective
% value. 
Conflict=abs(RewardArs-EffortArs);
% Add new column to datTable
DatTable_c=DatTable;
DatTable_c.Conflict = Conflict;
lmec = fitlme(DatTable_c,'PFC_theta~ Conflict + (Conflict|sbj_n_A)')
% Conflict here is going to be related to SV - high reward 

OV=RewardArs-EffortArs;
DatTable_co=DatTable_c;
DatTable_co.OV = OV;
lmeo = fitlme(DatTable_co,'PFC_beta~ OV + (OV|sbj_n_A)'); %OV p = 0.02

clc; close all;
DatTable_cod=DatTable_co;
DatTable_cod.DiffA = DiffA;
lmdiff = fitlme(DatTable_cod,'PFC_theta~ DiffA + (DiffA|ptnumA)')
