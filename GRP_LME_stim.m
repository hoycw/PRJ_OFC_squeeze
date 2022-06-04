%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
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
plot_time_lim = [-2 2];
stim_time_lim = [0 2];
plot_freq_lim = [2 30];

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
% Make a table

% Frequency & Time %
frt1=4;
frt2=7;

frb1=13;
frb2=30;
t1=0.5;
t2=1.5;

% Initialise variables
ptnumA=[];
datBmA=[];
datTmA=[];

StakeA=[];
RewardA=[];
EffortA=[];
EFFA=[];
DecisionA=[];
SubjValA=[];
DiffA=[];

RewardAs=[];
EffortAs=[];
EffortOAs=[];
DecisionAs=[];
SubjValAs=[];

for s=1:4
    
    PT=DataStr(s).data;
    
    frq=PT.TFbl.freq;
    tme=PT.TFbl.time;
    pow=PT.TFbl.powspctrm;
   
    bevar=PT.exp;
    
    % Run it on beta specific peaks %
    % Change
    btfr=[10,17,13,12];
    % Change with first one using peak %
%     btfr=[17,17,13,12];
    % Peak with 3rd one using change %
%     btfr=[17,22,13,13];
    
    frb1=btfr(s)-2;
    frb2=btfr(s)+2;
    
    % find index of frequency
    [mn,ixb1]=min((abs(frq-betapk_lim(s,1))));
    [mn,ixb2]=min((abs(frq-betapk_lim(s,2))));
    [mn,ixth1]=min((abs(frq-theta_lim(1))));
    [mn,ixth2]=min((abs(frq-theta_lim(2))));
    % find index of time
    [mn,ixt1]=min((abs(tme-t1)));
    [mn,ixt2]=min((abs(tme-t2)));
    % average beta across time, then frequencies
    datBext=squeeze(pow(:,:,ixb1:ixb2,ixt1:ixt2));
    datBm=mean(datBext,4);
    datBm=mean(datBm,3);
    % average theta across time, then frequencies
    datText=squeeze(pow(:,:,ixth1:ixth2,ixt1:ixt2));
    datTm=mean(datText,4);
    datTm=mean(datTm,3);
    
    % Concatenate SBJ, beta, theta values
    trlnum=size(datTm,1);
    
    ptnumber=num2str(ones(trlnum,1).*s);
    ptnumA=[ptnumA;ptnumber];
    
    datBmA=[datBmA;datBm];
    datTmA=[datTmA;datTm];
    
    % Extract the behavioral variables %
    RewardA=[RewardA; bevar.stake];
    EffortA=[EffortA; bevar.effort];
    EFFA=[EFFA; bevar.EFFs];
    DecisionA=[DecisionA; bevar.decision];
    SubjValA=[SubjValA; bevar.SV];
    DiffA=[DiffA;abs(bevar.SV-median(bevar.SV))];
    
    % Reward
    RewTmp=bevar.stake;
    RewTmp=smooth(RewTmp,2); % why isn't this doing anything?
    RewTmp(end-1:end)=[]; RewTmp=[nan(2,1); RewTmp]; % why shift by 2 instead of 1?
    RewardAs=[RewardAs; RewTmp];
    
    % Effort
    EffTmp=bevar.effort;
    EffTmp=smooth(EffTmp,1);
    EffTmp(end)=[]; EffTmp=[nan(1,1); EffTmp];
    EffortAs=[EffortAs; EffTmp];
    
    % Obj. Effort
    EffOTmp=bevar.EFFs;
    EffOTmp=smooth(EffOTmp,1);
    EffOTmp(end)=[]; EffOTmp=[nan(1,1); EffOTmp];
    EffortOAs=[EffortOAs; EffOTmp];
    
    % DecisionA
    DecTmp=bevar.decision;
    DecTmp=smooth(DecTmp,1);
    DecTmp(end)=[]; DecTmp=[nan(1,1); DecTmp];
    DecisionAs=[DecisionAs; DecTmp];
    
    % Subj Val.
    SVTmp=bevar.SV;
    SVTmp=smooth(SVTmp,1);
    SVTmp(end)=[]; SVTmp=[nan(1,1); SVTmp];
    SubjValAs=[SubjValAs; SVTmp];
    
end

OFCTheta=datTmA(:,2);
OFCBeta=datBmA(:,2);
LFPTheta=datTmA(:,1);
LFPBeta=datBmA(:,1);

DatTable = table(ptnumA,OFCTheta,OFCBeta,LFPTheta,LFPBeta,RewardA, EffortA, DecisionA, SubjValA, EFFA);
DatTable2 = table(ptnumA,OFCTheta,OFCBeta,LFPTheta,LFPBeta,RewardAs, EffortAs, DecisionAs, SubjValAs, EffortOAs);

DatTable2tmp=DatTable2;
DatTable2tmp.ptnumA=[];
DatTable2tmp.OFCTheta=[];
DatTable2tmp.OFCBeta=[];
DatTable2tmp.LFPTheta=[];
DatTable2tmp.LFPBeta=[];

DatTable3 = [DatTable, DatTable2tmp];


%% LME Modelliing %%
% OFC:
% Current trial
lme = fitlme(DatTable,'OFCBeta~ SubjValA + (1|ptnumA)')     %SV p = 0.01
lme = fitlme(DatTable,'OFCTheta~ SubjValA + (1|ptnumA)')    %SV p = 0.4

return;
% Previous trial
lme = fitlme(DatTable2,'OFCTheta~ SubjValAs + (1|ptnumA)')  %SV p = 0.04
lme = fitlme(DatTable2,'OFCBeta~ SubjValAs + (1|ptnumA)')   %
lme = fitlme(DatTable2,'LFPTheta~ SubjValAs + (1|ptnumA)')
lme = fitlme(DatTable2,'LFPBeta~ SubjValAs + (1|ptnumA)')


lme = fitlme(DatTable,'LFPBeta~ SubjValA + (1|ptnumA)')
lme = fitlme(DatTable,'LFPTheta~ SubjValA + (1|ptnumA)')


lme = fitlme(DatTable2,'OFCTheta~ SubjValAs + (1|ptnumA)')  %SV p = 0.04


% predict decision from neural data
glme = fitglme(DatTable3, 'DecisionA~ OFCBeta + OFCTheta + (1|ptnumA)','Link','log')    % beta p = 0.017
glme = fitglme(DatTable3, 'DecisionA~ LFPBeta + LFPTheta + (1|ptnumA)','Link','log')


return;

lme = fitlme(DatTable,'OFCTheta~ SubjValA + (SubjValA|ptnumA)')
lme = fitlme(DatTable,'OFCBeta~ SubjValA + (1|ptnumA)')
lme = fitlme(DatTable,'OFCBeta~ SubjValA*DecisionA + (SubjValA|ptnumA) + (DecisionA|ptnumA)')
lme = fitlme(DatTable2,'OFCTheta~ SubjValAs*DecisionAs + (SubjValAs|ptnumA) + (DecisionAs|ptnumA)')
lme = fitlme(DatTable2,'OFCTheta~ SubjValAs + (SubjValAs|ptnumA)')
lme = fitlme(DatTable2,'OFCTheta~ SubjValAs*DecisionAs + (1|ptnumA) + (1|ptnumA)')
lme = fitlme(DatTable2,'OFCTheta~ SubjValAs + (SubjValAs|ptnumA)')
lme = fitlme(DatTable2,'OFCTheta~ DecisionAs + (DecisionAs|ptnumA)')


% lme = fitlme(DatTable2,'OFCTheta~ DecisionAs + (DecisionAs|ptnumA)')
% 


return;

% Objective value%
% ObjVal=RewardA-EffortA

% lme = fitlme(DatTable,'OFCBeta~ SubjValA + (1|ptnumA)')
% lme = fitlme(DatTable,'OFCBeta~ SubjValA + (1|ptnumA) + (SubjValA-1|ptnumA)')
% figure();
% plotResiduals(lme,'fitted')
% find(residuals(lme) > 1.5)

lme2 = fitlme(DatTable,'OFCBeta ~ SubjValA ')
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
lmec = fitlme(DatTable_c,'OFCTheta~ Conflict + (Conflict|ptnumA)')
% Conflict here is going to be related to SV - high reward 

OV=RewardArs-EffortArs;
DatTable_co=DatTable_c;
DatTable_co.OV = OV;
lmeo = fitlme(DatTable_co,'OFCBeta~ OV + (OV|ptnumA)'); %OV p = 0.02

clc; close all;
DatTable_cod=DatTable_co;
DatTable_cod.DiffA = DiffA;
lmdiff = fitlme(DatTable_cod,'OFCTheta~ DiffA + (DiffA|ptnumA)')
