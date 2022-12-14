addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
clear all
close all

%%
model_p_acc = 0;
off_color = 'b';
on_color  = 'r';
scat_sz = 40;
font_sz = 18;

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
save_fig = 1;
fig_ftype = 'png';

%% Prepare the regression %
% This loads a 200x7 matrix derived from Behavioral_Summary.xlsx
% CleanReg.m- loads Reg.mat from xlsx, converts effort and stake to true
%   numbers instead of indices
% matrix columns: effort, stake, decision, block number (OFF: 1,4,6,7; ON: 2,3,5,8),
%   total trial, number, within block trial number, ON/OFF (1/0)
stim_bhv_fname = '/Users/colinhoy/Code/PRJ_OFC_squeeze/data/PFC05/Session2_withStimulation/Reg2.mat';
load(stim_bhv_fname);

data = Reg2;
col_names = {'Effort','Stake','Decision','Block','TrialCum','BlkCum','Stim'};

%% BEHAVIORAL MODELLING %%
effort_ix = strcmp(col_names,'Effort');
stake_ix  = strcmp(col_names,'Stake');
dec_ix    = strcmp(col_names,'Decision');
stim_ix   = strcmp(col_names,'Stim');
% Fit all behavior. Minimise the difference between the probability and the decision
decisionfun=@(p) norm( (exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2))))) - data(:,dec_ix));
[par_all, fit]=fminsearch(decisionfun, [1,1]);

SV_fn_all    = @(k) data(:,stake_ix)-(k*(data(:,effort_ix)).^2);
EFF_fn_all   = @(k) (k*(data(:,effort_ix)).^2);
SV_all   = SV_fn_all(par_all(2));
EFFs_all = EFF_fn_all(par_all(2));
p_accept = (exp(par_all(1)*(data(:,stake_ix)-(par_all(2)*(data(:,effort_ix)).^2))) ./...
    (exp(par_all(1)) + exp(par_all(1)*(data(:,stake_ix)-(par_all(2)*(data(:,effort_ix)).^2)))));

% Fit OFF behavior
off_idx = data(:,stim_ix)==0;
decisionfun=@(p) norm( (exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2))))) - data(off_idx,dec_ix));
[par_off, fit]=fminsearch(decisionfun, [1,1]);

p_accept_off = (exp(par_off(1)*(data(off_idx,stake_ix)-(par_off(2)*(data(off_idx,effort_ix)).^2))) ./...
    (exp(par_off(1)) + exp(par_off(1)*(data(off_idx,stake_ix)-(par_off(2)*(data(off_idx,effort_ix)).^2)))));
SV_fn_off    = @(k) data(off_idx,stake_ix)-(k*(data(off_idx,effort_ix)).^2);
SV_off   = SV_fn_off(par_off(2));
[~,SV_off_sort_idx] = sort(SV_off);

% Fit ON behavior
on_idx = data(:,stim_ix)==1;
decisionfun=@(p) norm( (exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2))))) - data(on_idx,dec_ix));
[par_on, fit]=fminsearch(decisionfun, [1,1]);

p_accept_on = (exp(par_on(1)*(data(on_idx,stake_ix)-(par_on(2)*(data(on_idx,effort_ix)).^2))) ./...
    (exp(par_on(1)) + exp(par_on(1)*(data(on_idx,stake_ix)-(par_on(2)*(data(on_idx,effort_ix)).^2)))));
SV_fn_on    = @(k) data(on_idx,stake_ix)-(k*(data(on_idx,effort_ix)).^2);
SV_on   = SV_fn_on(par_on(2));
[~,SV_on_sort_idx] = sort(SV_on);

%% Fit a Generalised linear model and look for an interaction term %%
zReg = [zscore(data(:,1)) zscore(data(:,2)) data(:,3:7)];
if model_p_acc
    zReg = [zReg zscore([p_accept_off; p_accept_on])];
    col_names = [col_names 'p_accept'];
end
% T = array2table(data, 'VariableNames', col_names);
zT = array2table(zReg, 'VariableNames', col_names);


% GLMs
% modelspec = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim';
% mdl = fitglm(T,modelspec,'Distribution','binomial')
% zmdl = fitglm(zT,modelspec,'Distribution','binomial')

% GLMMs
modelspec = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
% mdl = fitglme(T,modelspec,'Distribution','binomial')
zmdl = fitglme(zT,modelspec,'Distribution','binomial')

% GLMMs
if model_p_acc
    modelspec = 'p_accept ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
    zmdl_pacc = fitglme(zT,modelspec)
end

% mdl2 = stepwiselm(T,'interactions')
b1=[];
beta=zmdl.Coefficients.Estimate

b1=beta(6:7)

% addpath('E:\Box Sync\Research\Resources & References\Software\Matlab\General_Code');

%% Average behavior by condition for plotting
% Mean Reward and Effort
efforts = unique(data(:,effort_ix));
stakes  = unique(data(:,stake_ix));
if length(efforts)~=5 || length(stakes)~=5; error('should be 5 conditions!'); end
p_dec_mn     = nan([length(stakes) length(efforts)]);
p_dec_mn_on  = nan([length(stakes) length(efforts)]);
p_dec_mn_off = nan([length(stakes) length(efforts)]);
p_acc_mn     = nan([length(stakes) length(efforts)]);
p_acc_mn_on  = nan([length(stakes) length(efforts)]);
p_acc_mn_off = nan([length(stakes) length(efforts)]);
dec_stake_off_mn = nan([length(stakes) 1]);
dec_stake_off_se = nan([length(stakes) 1]);
dec_stake_on_mn  = nan([length(stakes) 1]);
dec_stake_on_se  = nan([length(stakes) 1]);
pacc_stake_off_mn = nan([length(stakes) 1]);
pacc_stake_off_se = nan([length(stakes) 1]);
pacc_stake_on_mn  = nan([length(stakes) 1]);
pacc_stake_on_se  = nan([length(stakes) 1]);
for stk_i = 1:5
    % OFF stim: Average decision by stake
    cond_idx = data(:,stake_ix)==stakes(stk_i) & off_idx;
    dec_stake_off_mn(stk_i) = mean(data(cond_idx,dec_ix));
    dec_stake_off_se(stk_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % OFF stim: Average probability of acceptance by stake
    trl_idx = cond_idx(off_idx);
    pacc_stake_off_mn(stk_i) = mean(p_accept_off(trl_idx));
    pacc_stake_off_se(stk_i) = std(p_accept_off(trl_idx))./sqrt(sum(trl_idx));
    
    % ON stim: Average decision by stake
    cond_idx = data(:,stake_ix)==stakes(stk_i) & on_idx;
    dec_stake_on_mn(stk_i) = mean(data(cond_idx,dec_ix));
    dec_stake_on_se(stk_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % ON stim: Average probability of acceptance by stake
    trl_idx = cond_idx(on_idx);
    pacc_stake_on_mn(stk_i) = mean(p_accept_on(trl_idx));
    pacc_stake_on_se(stk_i) = std(p_accept_on(trl_idx))./sqrt(sum(trl_idx));
    
    for eff_i = 1:5
        cond_idx = data(:,effort_ix)==efforts(eff_i) & data(:,stake_ix)==stakes(stk_i);
        % Average over all data
        p_dec_mn(stk_i,eff_i) = mean(data(cond_idx,dec_ix));
        p_acc_mn(stk_i,eff_i) = mean(p_accept(cond_idx));
        % Average over OFF data
        trl_idx = cond_idx(off_idx);
        p_dec_mn_off(stk_i,eff_i) = mean(data(cond_idx & off_idx,dec_ix));
        p_acc_mn_off(stk_i,eff_i) = mean(p_accept_off(trl_idx));
        % Average over ON data
        trl_idx = cond_idx(on_idx);
        p_dec_mn_on(stk_i,eff_i) = mean(data(cond_idx & on_idx,dec_ix));
        p_acc_mn_on(stk_i,eff_i) = mean(p_accept_on(trl_idx));
    end
end

% Prob accept ON - OFF
p_acc_mn_diff = p_acc_mn_on-p_acc_mn_off;
p_dec_mn_diff = p_dec_mn_on-p_dec_mn_off;

%% Line plots of Prob accept by reward and effort
fig_name = ['PFC05_stim_allTrials_dec_stake_ONOFF_line'];
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(stakes,dec_stake_off_mn,dec_stake_off_se,'Color',off_color,'linewidth',3);
errorbar(stakes,dec_stake_on_mn,dec_stake_on_se,'Color',on_color,'linewidth',3);
xlabel('Reward');
ylabel('% Accept');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Decision ~ Stimulation*Reward (p=' num2str(zmdl.Coefficients.pValue(7),'%.03f') ')']);
set(gca,'FontSize',18);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/decision_stake_line/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
%     fig_fname = [fig_dir fig_name '.fig'];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
end

fig_name = ['PFC05_stim_allTrials_p_acc_stake_ONOFF_line'];
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(stakes,pacc_stake_off_mn,pacc_stake_off_se,'Color',off_color,'linewidth',3);
errorbar(stakes,pacc_stake_on_mn,pacc_stake_on_se,'Color',on_color,'linewidth',3);
xlabel('Reward');
ylabel('Probability of Accept');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Effect of Stimulation*Reward on Prob. Accept']);
set(gca,'FontSize',18);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
%     fig_fname = [fig_dir fig_name '.fig'];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
end

%% Line plot for SV vs. decision function
fig_name = ['PFC05_stim_dec_fn_ONOFF_line'];
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
off_line = line(SV_off(SV_off_sort_idx),p_accept_off(SV_off_sort_idx),'Color',off_color,'LineWidth',3);
scatter(SV_off(SV_off_sort_idx),p_accept_off(SV_off_sort_idx),scat_sz,off_color);
on_line = line(SV_on(SV_on_sort_idx),p_accept_on(SV_on_sort_idx),'Color',on_color,'LineWidth',3);
scatter(SV_on(SV_on_sort_idx),p_accept_on(SV_on_sort_idx),scat_sz,on_color);
%     line(mdl_x,mdl_y,'Color','k');
line([0 0],[0 1],'Color','k','LineStyle','--');
xlabel('Subjective Value'); ylabel('Probabilty Accept');
title('Subjective Value Decision Function');
legend([off_line,on_line],{'OFF stimulation','ON stimulation'},'Location','best');
set(gca,'FontSize',font_sz);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/decision_fn_line/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
%     fig_fname = [fig_dir fig_name '.fig'];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
end

%% 3-D plot of decision as function of reward and effort
% Total behavior
fig_name = ['PFC05_stim_allTrials_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_dec_mn);
zlabel('% Accept');
title('PFC05 Decision Landscape: All trials');
caxis([0 1]);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/p_decision_bar3/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON-OFF behavior
fig_name = ['PFC05_stim_ON-OFF_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_dec_mn_diff);
zlabel('% Accept');
title('PFC05 Decision Landscape: ON-OFF');
rb_cmap = redblue();
colormap(rb_cmap);
caxis([-max(abs(p_dec_mn_diff(:))) max(abs(p_dec_mn_diff(:)))]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% OFF stim
fig_name = ['PFC05_stim_OFF_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_dec_mn_off);
zlabel('% Accept');
title('PFC05 Decision Landscape: OFF Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON stim
fig_name = ['PFC05_stim_ON_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_dec_mn_on);
zlabel('% Accept');
title('PFC05 Decision Landscape: ON Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% 3-D plot of acceptance as function of reward and effort
% Total behavior
fig_name = ['PFC05_stim_allTrials_p_acc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_acc_mn);

title('PFC05 Decision Model Landscape: All trials');
caxis([0 1]);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/p_accept_bar3/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON-OFF behavior
fig_name = ['PFC05_stim_ON-OFF_p_acc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_acc_mn_diff);
title('PFC05 Decision Model Landscape: ON-OFF');
rb_cmap = redblue();
colormap(rb_cmap);
caxis([-max(abs(p_acc_mn_diff(:))) max(abs(p_acc_mn_diff(:)))]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% OFF stim
fig_name = ['PFC05_stim_OFF_p_acc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_acc_mn_off);
title('PFC05 Decision Model Landscape: OFF Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON stim
fig_name = ['PFC05_stim_ON_p_acc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,p_acc_mn_on);
title('PFC05 Decision Model Landscape: ON Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% define function for 3D bar plot
function [bars] = fn_3d_decision_landscape(stakes,efforts,zvals)

bars = bar3(zvals);
ylabel('Reward');
yticklabels(stakes);
xlabel('Effort (% Max)');
xticklabels(efforts*100);
zlabel('Probability Accept');
set(gca,'FontSize',font_sz);
view([135 30]);

% Shade bar according to z-axis
for b = 1:length(bars)
    zdata = bars(b).ZData;
    bars(b).CData = zdata;
    bars(b).FaceColor = 'interp';
end
%     idxrep = repmat(1:size(p_acc_mn_stim,1),6,1);  %magic 6
%     for col_ix = 1:size(p_acc_mn_stim,2)
%         Z = repmat(p_acc_mn_stim(idxrep(:),col_ix,p),1,4);   %magic 4
%         set(bars(col_ix), 'CData',Z);
%     end
end


%% ========================================================================
%           OLD SIMON CODE
%  ========================================================================
% %% Contrast-based GLM
% % Organise Dependent and Independent Variables %
% % Reg
% % 1=force 2=reward 3 = decision (1= accept; 0 = reject)7 = stim 
% 
% Y = data(:,3);
% % X = Reg(:,[1:2 5 6 7]);
% X = data(:,[1:2 7]);
% 
% Xes =data(:,1).*data(:,7);
% Xrs =data(:,2).*data(:,7);
% % Add a line of ones %
% 
% SV=Xrs-Xes;
% 
% Xn=[X, Xes, Xrs]
% % 
% 
% [B,dev,stats] = glmfit(Xn,Y,'binomial', 'link', 'logit');
% B
% stats.p
% 
% %%
% figure;
% set(gcf, 'Position',  [100, 100, 180, 180]);
% b=bar(b1,'FaceColor', rgb('LightSeaGreen'))
% box off
% 
% % Aquamarine LightSeaGreen Teal
% 
% %% Get absolute behavior
% 
% wh1=find(data(:,7)==1);
% wh0=find(data(:,7)==0);
% 
% Ym1(1)=mean(Y(wh1))
% Ym1(2)=mean(Y(wh0))
% 
% Ym1=Ym1*100;
% figure;
% set(gcf, 'Position',  [100, 100, 180, 180]);
% b=bar(Ym1,'FaceColor','flat'); ylim([70 80]);
% b.CData(2,:) = 'r';%rgb('Red');
% box off