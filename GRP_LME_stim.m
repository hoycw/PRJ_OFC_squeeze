%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
SBJ_pfc_rois  = {'FPC', 'OFC', 'OFC', 'FPC'};
% bg_roi   = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_colors = distinguishable_colors(length(SBJs));

% Analysis parameters:
theta_lim  = [4 7];
beta_lim   = [13 30];    
sbj_beta_pk = [10,17,13,12]; % PFC03, PFC04, PFC05, PFC01
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];
betapk_bw = 4;
betapk_lim = nan(length(SBJs),2);
for s = 1:length(sbj_beta_pk)
    betapk_lim(s,:) = [sbj_beta_pk(s)-betapk_bw/2 sbj_beta_pk(s)+betapk_bw/2];
end

% Plotting parameters
sem_alpha  = 0.5;
time_lim = [0.5 1.5];
% stim_time_lim = [0 2];
% plot_freq_lim = [2 30];

save_fig = 0;

%% Load data
% clc
% close all
% clear all

% restoredefaultpath
% addpath('C:\BackupRepo\fieldtrip-20190705');
% addpath('C:\Users\slittle\Box\Research\Resources & References\Software\Matlab\General_Code');
% ft_defaults

DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';

[numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);
SyncDetails = numbers;

% Load data files
DataStr=struct;
for s=1:size(FileDetails,1)-1
    
    load(strcat(DataStorage2,'/',FileDetails{s,1},'Stimulus_Locked.mat'));
    DataStr(s).data=AllData;
end

%% Convert into format that is suitable for LME modelling %%
% Initialize channel, frequency, and time variables
frqs = DataStr(1).data.TFbl.freq;   % 2:80 Hz
tmes = DataStr(1).data.TFbl.time;   % -3:3 s
ch_lab = DataStr(1).data.TFbl.label;
ofc_ch_ix = find(strcmp(DataStr(1).data.TFbl.label,'OFC'));
lfp_ch_ix = find(strcmp(DataStr(1).data.TFbl.label,'LFP'));

% find index of time
[~,tm_ix(1)] = min(abs(tmes-time_lim(1)));
[~,tm_ix(2)]    = min(abs(tmes-time_lim(2)));

% find frequency indices
for i = 1:2
    [~,th_ix(i)] = min(abs(frqs-theta_lim(i)));
    [~,b_ix(i)] = min(abs(frqs-beta_lim(i)));
    for s = 1:length(SBJs)
        [~,bpk_ix(s,i)] = min(abs(frqs-betapk_lim(s,i)));
    end
end


% Make a table

% Initialise variables (A = current trial, As = previous trial)
sbj_n_A = [];
pfc_roi = [];
trl_n_A = [];
bpk_A   = [];
th_A    = [];

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
    
    pow=DataStr(s).data.TFbl.powspctrm;
   
    bhv=DataStr(s).data.exp;
    
    % Run it on beta specific peaks %
    
    % average beta across time, then frequencies
    bpk_mn = mean(pow(:,:,bpk_ix(1):bpk_ix(2),tm_ix(1):tm_ix(2)),4);
    bpk_mn = mean(bpk_mn,3);
    % average theta across time, then frequencies
    th_mn = mean(pow(:,:,th_ix(1):th_ix(2),tm_ix(1):tm_ix(2)),4);
    th_mn = mean(th_mn,3);
    
    % Concatenate SBJ, beta, theta values
    trl_n = size(th_mn,1);
    trl_n_A = [trl_n_A; [1:trl_n]'];
    sbj_n_A = [sbj_n_A; num2str(ones(trl_n,1).*s)];
    if strcmp(SBJ_pfc_rois{s},'OFC'); roi_ix = 1; else; roi_ix = 2; end
    pfc_roi = [pfc_roi; num2str(ones(trl_n,1).*roi_ix)];
    
    bpk_A   = [bpk_A; bpk_mn];
    th_A    = [th_A; th_mn];
    
    % Extract the behavioral variables %
    reward_A   = [reward_A; bhv.stake];
    effort_A   = [effort_A; bhv.effort];
    EFF_A      = [EFF_A; bhv.EFFs];
    decision_A = [decision_A; bhv.decision];
    SV_A       = [SV_A; bhv.SV];
    SV_md_A    = [SV_md_A; abs(bhv.SV-median(bhv.SV))];
    
    % Reward
    rew_prv = smooth(bhv.stake,2); % why isn't this doing anything?
    rew_prv(end-1:end)=[]; rew_prv=[nan(2,1); rew_prv]; % why shift by 2 instead of 1?
    reward_As = [reward_As; rew_prv];
    
    % Effort
    effort_prv = smooth(bhv.effort,1);
    effort_prv(end)=[]; effort_prv = [nan(1,1); effort_prv];
    effort_As = [effort_As; effort_prv];
    
    % Obj. Effort
    effortO_prv = smooth(bhv.EFFs,1);
    effortO_prv(end)=[]; effortO_prv = [nan(1,1); effortO_prv];
    effortO_As = [effortO_As; effortO_prv];
    
    % decision_A
    dec_prv = smooth(bhv.decision,1);
    dec_prv(end) = []; dec_prv = [nan(1,1); dec_prv];
    decision_As = [decision_As; dec_prv];
    
    % Subj Val.
    SV_prv = smooth(bhv.SV,1);
    SV_prv(end) = []; SV_prv = [nan(1,1); SV_prv];
    SV_As = [SV_As; SV_prv];
end

OFC_theta = th_A(:,ofc_ch_ix);
OFC_beta  = bpk_A(:,ofc_ch_ix);
LFP_theta = th_A(:,lfp_ch_ix);
LFP_beta  = bpk_A(:,lfp_ch_ix);

table_A  = table(trl_n_A, sbj_n_A, pfc_roi, OFC_theta, OFC_beta, LFP_theta, LFP_beta,...
                 reward_A, effort_A, decision_A, SV_A, EFF_A);
table_As = table(trl_n_A, sbj_n_A, pfc_roi, OFC_theta, OFC_beta, LFP_theta, LFP_beta,...
                 reward_As, effort_As, decision_As, SV_As, effortO_As);

table_As_only = table(reward_As, effort_As, decision_As, SV_As, effortO_As);

DatTable2tmp=table_As;
DatTable2tmp.sbj_n_A=[];
DatTable2tmp.trl_n_A=[];
DatTable2tmp.pfc_roi=[];
DatTable2tmp.OFC_theta=[];
DatTable2tmp.OFC_beta=[];
DatTable2tmp.LFP_theta=[];
DatTable2tmp.LFP_beta=[];

DatTable3 = [table_A, DatTable2tmp];
table_all = [table_A, table_As_only];
% These are different by trials 1 and 2 for each patient, unclear why (they
% look identical...)

%% LME Modelliing %%
% OFC:
% Current trial
lme = fitlme(table_A,'OFC_beta~ SV_A + (1|sbj_n_A)')     %SV p = 0.01
lme = fitlme(table_A,'OFC_theta~ SV_A + (1|sbj_n_A)')    %SV p = 0.4

lme = fitlme(table_A,'OFC_beta~ SV_A*pfc_roi + (1|sbj_n_A)')     %SV p = 0.0007; pfc_roi p = 0.0048; * p = 0.1
% lme = fitlme(table_A,'OFC_beta~ SV_A + (SV_A|pfc_roi) + (1|sbj_n_A)')     %SV p = 0.048
lme = fitlme(table_A,'OFC_theta~ SV_A*pfc_roi + (1|sbj_n_A)')    %none

return;
% Previous trial
lme = fitlme(table_As,'OFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04
lme = fitlme(table_As,'OFC_beta~ SV_As + (1|sbj_n_A)')   %
lme = fitlme(table_As,'LFP_theta~ SV_As + (1|sbj_n_A)')
lme = fitlme(table_As,'LFP_beta~ SV_As + (1|sbj_n_A)')

lme = fitlme(table_As,'OFC_theta~ SV_As*pfc_roi + (1|sbj_n_A)')  %SV p = 0.052
lme = fitlme(table_As,'OFC_theta~ SV_As + (SV_As|pfc_roi) + (1|sbj_n_A)')  %SV p = 0.038

lme = fitlme(table_A,'LFP_beta~ SV_A + (1|sbj_n_A)')
lme = fitlme(table_A,'LFP_theta~ SV_A + (1|sbj_n_A)')


lme = fitlme(table_As,'OFC_theta~ SV_As + (1|sbj_n_A)')  %SV p = 0.04


% predict decision from neural data
glme = fitglme(DatTable3, 'decision_A~ OFC_beta + OFC_theta + (1|sbj_n_A)','Link','log')    % beta p = 0.017
glme = fitglme(DatTable3, 'decision_A~ LFP_beta + LFP_theta + (1|sbj_n_A)','Link','log')


return;

lme = fitlme(DatTable,'OFC_theta~ SV_A + (SV_A|sbj_n_A)')
lme = fitlme(DatTable,'OFC_beta~ SV_A + (1|sbj_n_A)')
lme = fitlme(DatTable,'OFC_beta~ SV_A*decision_A + (SV_A|sbj_n_A) + (decision_A|sbj_n_A)')
lme = fitlme(table_As,'OFC_theta~ SV_As*decision_As + (SV_As|sbj_n_A) + (decision_As|sbj_n_A)')
lme = fitlme(table_As,'OFC_theta~ SV_As + (SV_As|sbj_n_A)')
lme = fitlme(table_As,'OFC_theta~ SV_As*decision_As + (1|sbj_n_A) + (1|sbj_n_A)')
lme = fitlme(table_As,'OFC_theta~ SV_As + (SV_As|sbj_n_A)')
lme = fitlme(table_As,'OFC_theta~ decision_As + (decision_As|sbj_n_A)')


% lme = fitlme(table_As,'OFC_theta~ decision_As + (decision_As|sbj_n_A)')
% 


return;

% Objective value%
% ObjVal=RewardA-EffortA

% lme = fitlme(DatTable,'OFC_beta~ SV_A + (1|sbj_n_A)')
% lme = fitlme(DatTable,'OFC_beta~ SV_A + (1|sbj_n_A) + (SV_A-1|sbj_n_A)')
% figure();
% plotResiduals(lme,'fitted')
% find(residuals(lme) > 1.5)

lme2 = fitlme(DatTable,'OFC_beta ~ SV_A ')
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
lmec = fitlme(DatTable_c,'OFC_theta~ Conflict + (Conflict|sbj_n_A)')
% Conflict here is going to be related to SV - high reward 

OV=RewardArs-EffortArs;
DatTable_co=DatTable_c;
DatTable_co.OV = OV;
lmeo = fitlme(DatTable_co,'OFC_beta~ OV + (OV|sbj_n_A)'); %OV p = 0.02

clc; close all;
DatTable_cod=DatTable_co;
DatTable_cod.DiffA = DiffA;
lmdiff = fitlme(DatTable_cod,'OFC_theta~ DiffA + (DiffA|ptnumA)')
