clear all
close all

%% Prepare the regression %

stim_bhv_fname = '/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/PFC05/Session2_withStimulation/Reg2.mat';
load(stim_bhv_fname);

Reg=Reg2;
%  

%%
% Organise Dependent and Independent Variables %
% Reg
% 1=force 2=reward 3 = decision (1= accept; 0 = reject)7 = stim 

Y = Reg(:,3);
% X = Reg(:,[1:2 5 6 7]);
X = Reg(:,[1:2 7]);

Xes =Reg(:,1).*Reg(:,7);
Xrs =Reg(:,2).*Reg(:,7);
% Add a line of ones %

SV=Xrs-Xes;

Xn=[X, Xes, Xrs]
% 

[B,dev,stats] = glmfit(Xn,Y,'binomial', 'link', 'logit');
B
stats.p


%% Fit a Generalised linear model and look for an interaction term %%
zReg = [zscore(Reg(:,1)) zscore(Reg(:,2)) Reg(:,3:7)];
T = array2table(Reg, 'VariableNames', {'Effort','Stake','Decision','Block','TrialCum','BlkCum','Stim'})
zT = array2table(zReg, 'VariableNames', {'Effort','Stake','Decision','Block','TrialCum','BlkCum','Stim'})

modelspec = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim';

% GLMs
mdl = fitglm(T,modelspec,'Distribution','binomial')
zmdl = fitglm(zT,modelspec,'Distribution','binomial')

% GLMMs
modelspec = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
mdl = fitglme(T,modelspec,'Distribution','binomial')
zmdl = fitglme(zT,modelspec,'Distribution','binomial')


% mdl2 = stepwiselm(T,'interactions')
b1=[];
beta=mdl.Coefficients.Estimate

b1=beta(6:7)

addpath('E:\Box Sync\Research\Resources & References\Software\Matlab\General_Code');

%%
figure;
set(gcf, 'Position',  [100, 100, 180, 180]);
b=bar(b1,'FaceColor', rgb('LightSeaGreen'))
box off

% Aquamarine LightSeaGreen Teal

%% Get absolute behavior

wh1=find(Reg(:,7)==1);
wh0=find(Reg(:,7)==0);

Ym1(1)=mean(Y(wh1))
Ym1(2)=mean(Y(wh0))

Ym1=Ym1*100;
figure;
set(gcf, 'Position',  [100, 100, 180, 180]);
b=bar(Ym1,'FaceColor','flat'); ylim([70 80]);
b.CData(2,:) = rgb('Red');
box off