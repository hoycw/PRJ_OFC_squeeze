addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
clear all
close all

%%
model_pacc = 1;
off_color = 'b';
on_color  = 'r';
scat_sz = 40;
font_sz = 24;

stim_cond = {'OFF','ON'};

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
save_fig = 1;
fig_ftype = 'svg';%'png';

%% Prepare the regression %
% This loads a 200x7 matrix derived from Behavioral_Summary.xlsx
% CleanReg.m- loads Reg.mat from xlsx, converts effort and stake to true
%   numbers instead of indices
% matrix columns: effort, stake, decision, block number (OFF: 1,4,6,7; ON: 2,3,5,8),
%   total trial, number, within block trial number, ON/OFF (1/0)
stim_bhv_fname = '/Users/colinhoy/Code/PRJ_OFC_squeeze/data/PFC05/Session2_withStimulation/Reg2.mat';
load(stim_bhv_fname);

data = Reg2;
col_names = {'Effort','Stake','Decision','Block','TrialCum','BlkCum','Stim','p_accept'};

%% BEHAVIORAL MODELLING %%
effort_ix = strcmp(col_names,'Effort');
stake_ix  = strcmp(col_names,'Stake');
dec_ix    = strcmp(col_names,'Decision');
stim_ix   = strcmp(col_names,'Stim');
blk_ix    = strcmp(col_names,'Block');
% Fit all behavior. Minimise the difference between the probability and the decision
decisionfun=@(p) norm( (exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2))))) - data(:,dec_ix));
[par_all, fit_all]=fminsearch(decisionfun, [1,1]);

SV_fn_all    = @(k) data(:,stake_ix)-(k*(data(:,effort_ix)).^2);
EFF_fn_all   = @(k) (k*(data(:,effort_ix)).^2);
SV_all   = SV_fn_all(par_all(2));
EFFs_all = EFF_fn_all(par_all(2));
pacc = (exp(par_all(1)*(data(:,stake_ix)-(par_all(2)*(data(:,effort_ix)).^2))) ./...
    (exp(par_all(1)) + exp(par_all(1)*(data(:,stake_ix)-(par_all(2)*(data(:,effort_ix)).^2)))));

% Fit OFF behavior
off_idx = data(:,stim_ix)==0;
decisionfun=@(p) norm( (exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2))))) - data(off_idx,dec_ix));
[par_off, fit_off]=fminsearch(decisionfun, [1,1]);

pacc_off = (exp(par_off(1)*(data(off_idx,stake_ix)-(par_off(2)*(data(off_idx,effort_ix)).^2))) ./...
    (exp(par_off(1)) + exp(par_off(1)*(data(off_idx,stake_ix)-(par_off(2)*(data(off_idx,effort_ix)).^2)))));
SV_fn_off    = @(k) data(off_idx,stake_ix)-(k*(data(off_idx,effort_ix)).^2);
SV_off   = SV_fn_off(par_off(2));
[~,SV_off_sort_idx] = sort(SV_off);

% Fit ON behavior
on_idx = data(:,stim_ix)==1;
decisionfun=@(p) norm( (exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2))))) - data(on_idx,dec_ix));
[par_on, fit_on]=fminsearch(decisionfun, [1,1]);

pacc_on = (exp(par_on(1)*(data(on_idx,stake_ix)-(par_on(2)*(data(on_idx,effort_ix)).^2))) ./...
    (exp(par_on(1)) + exp(par_on(1)*(data(on_idx,stake_ix)-(par_on(2)*(data(on_idx,effort_ix)).^2)))));
SV_fn_on    = @(k) data(on_idx,stake_ix)-(k*(data(on_idx,effort_ix)).^2);
SV_on   = SV_fn_on(par_on(2));
[~,SV_on_sort_idx] = sort(SV_on);

%% Add model with indifference point
% from Klein-Flugge 2015 PLoSCB:
%   p(Accept) = 1/(1 + exp(-B*(SV - a)) );
%   where a is a free parameter for the indifference point (sigmoid crosses 0.5)
dec_fun_ind=@(p) norm( (exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2) - p(3))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(:,stake_ix)-(p(2)*(data(:,effort_ix)).^2) - p(3))))) - data(:,dec_ix));
[par_ind, fit_ind]=fminsearch(dec_fun_ind, [1,1,1]);

% Run separately for OFF vs ON stim
dec_fun_ind_off=@(p) norm( (exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2) - p(3))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(off_idx,stake_ix)-(p(2)*(data(off_idx,effort_ix)).^2) - p(3))))) - data(off_idx,dec_ix));
[par_ind_off, fit_ind_off]=fminsearch(dec_fun_ind_off, [1,1,1]);

dec_fun_ind_on=@(p) norm( (exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2) - p(3))) ./ ...
    (exp(p(1)) + exp(p(1)*(data(on_idx,stake_ix)-(p(2)*(data(on_idx,effort_ix)).^2) - p(3))))) - data(on_idx,dec_ix));
[par_ind_on, fit_ind_on]=fminsearch(dec_fun_ind_on, [1,1,1]);

% Compute p(Accept)
pacc_ind_off = (exp(par_ind_off(1)*(data(off_idx,stake_ix)-(par_ind_off(2)*(data(off_idx,effort_ix)).^2) - par_ind_off(3))) ./ ...
    (exp(par_ind_off(1)) + exp(par_ind_off(1)*(data(off_idx,stake_ix)-(par_ind_off(2)*(data(off_idx,effort_ix)).^2) - par_ind_off(3)))));
SV_ind_off   = SV_fn_off(par_ind_off(2));
[~,SV_ind_off_sort_idx] = sort(SV_ind_off);
pacc_ind_on = (exp(par_ind_on(1)*(data(on_idx,stake_ix)-(par_ind_on(2)*(data(on_idx,effort_ix)).^2) - par_ind_on(3))) ./ ...
    (exp(par_ind_on(1)) + exp(par_ind_on(1)*(data(on_idx,stake_ix)-(par_ind_on(2)*(data(on_idx,effort_ix)).^2) - par_ind_on(3)))));
SV_ind_on   = SV_fn_on(par_ind_on(2));
[~,SV_ind_on_sort_idx] = sort(SV_ind_on);
sv_vals  = -4.5:0.01:13;
pacc_vals_on_ind  = (exp(par_ind_on(1)*sv_vals)) ./ (exp(par_ind_on(1)) + exp(par_ind_on(1)*sv_vals));
pacc_vals_off_ind = (exp(par_ind_off(1)*sv_vals)) ./ (exp(par_ind_off(1)) + exp(par_ind_off(1)*sv_vals));

% Plot the outcomes
fig_name = 'PFC05_stim_dec_fn_ind_ONOFF_line';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
off_line = plot(SV_ind_off(SV_ind_off_sort_idx),pacc_ind_off(SV_ind_off_sort_idx),'Color',off_color,'LineWidth',2);
on_line  = plot(SV_ind_on(SV_ind_on_sort_idx),pacc_ind_on(SV_ind_on_sort_idx),'Color',on_color,'LineWidth',2);
% off_line = plot(sv_vals,pacc_vals_off_ind,'Color',off_color,'LineWidth',2);
% on_line  = plot(sv_vals,pacc_vals_on_ind,'Color',on_color,'LineWidth',2);
legend([off_line, on_line],stim_cond);
line([-5 15],[0.5 0.5],'Color','k','LineStyle','--');
xlim([-5 15]);
xlabel('Subjective Value'); ylabel('Model Probabilty Accept');
title('Subjective Value Decision Function');
legend([off_line,on_line],strcat(stim_cond,' Stimulation'),'Location','best');
set(gca,'FontSize',font_sz);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/decision_fn_line_indif_pt/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% pacc_on_ind = (exp(par_ind_on(1)*par_ind_on(1)) ./ (exp(par_ind_on(1)) + exp(par_ind_on(1)*par_ind_on(3))));
% pacc_off_ind = (exp(par_ind_off(1)*par_ind_off(3)) ./ (exp(par_ind_off(1)) + exp(par_ind_off(1)*par_ind_off(3))));

%% Compile table
pacc_all = [pacc_off; pacc_on];
ease_all = abs(pacc_all - 0.5);
zReg = [zscore(data(:,effort_ix)) zscore(data(:,stake_ix)) data(:,3:7) zscore(pacc_all)];
% T = array2table(data, 'VariableNames', col_names);
tbl_z = array2table(zReg, 'VariableNames', col_names);

% Derive previous trial predictors
tbl_z.Stake_prv  = [nan; tbl_z.Stake(1:end-1)];
tbl_z.Effort_prv = [nan; tbl_z.Effort(1:end-1)];

% % Remove first trial of each block
% tbl_z(tbl_z.BlkCum==1,:) = array2table(nan);

%% Fit a Generalised linear model and look for an interaction term %%
% GLMMs
fit_method = 'Laplace';
full_mdl = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim';% + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
% lme = fitglme(T,modelspec,'Distribution','binomial','FitMethod',fit_method);)
full_lme = fitglme(tbl_z,full_mdl,'Distribution','binomial','FitMethod',fit_method);

% % Test if block helps
% full_mdl_blk = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block)';
% blk_lme = fitglme(tbl_z,full_mdl_blk,'Distribution','binomial','FitMethod',fit_method);
% blk_pval = compare(full_lme,blk_lme,'CheckNesting',true)
% % Nope, adding block is not better

% Test significance
preds = {'Stake','Effort','Stim','Effort:Stake','Stake:Stim','Effort:Stim'};
[dec_null_mdls, dec_comparisons] = fn_GLMM_run_null_LRT_fixed_effects(tbl_z,full_lme,preds,fit_method);
fprintf('Full Model: %s\n',full_lme.Formula);
for p_ix = 1:length(preds)
    fprintf('Reduced Model for %s: %s\n',preds{p_ix},dec_null_mdls{p_ix}.Formula);
    coef_ix = strcmp(full_lme.CoefficientNames,preds{p_ix});
    fprintf('%s: B = %.4f (p = %.4f)\n',preds{p_ix},full_lme.Coefficients.Estimate(coef_ix),dec_comparisons{p_ix}.pValue(2));
end
fprintf('\n\n');

% GLMMs
if model_pacc
    fit_method = 'ML';
    full_pacc_mdl = 'p_accept ~ Effort*Stake*Stim - Effort:Stake:Stim';% + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
    pacc_lme = fitlme(tbl_z,full_pacc_mdl,'FitMethod',fit_method);

    preds = {'Stake','Effort','Stim','Effort:Stake','Stake:Stim','Effort:Stim'};
    [pacc_null_mdls, pacc_comparisons] = fn_LMM_run_null_LRT_fixed_effects(tbl_z,pacc_lme,preds,fit_method);
    fprintf('Full Model: %s\n',pacc_lme.Formula);
    for p_ix = 1:length(preds)
        fprintf('Reduced Model for %s: %s\n',preds{p_ix},pacc_null_mdls{p_ix}.Formula);
        coef_ix = strcmp(pacc_lme.CoefficientNames,preds{p_ix});
        fprintf('%s: B = %.4f (p = %.4f)\n',preds{p_ix},pacc_lme.Coefficients.Estimate(coef_ix),pacc_comparisons{p_ix}.pValue(2));
    end
end

%% Try with previous trial predictors
% fit_method = 'Laplace';
% full_prv_mdl = 'Decision ~ Effort*Stake*Stim - Effort:Stake:Stim + Stake_prv + Effort_prv';% + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
% % lme = fitglme(T,modelspec,'Distribution','binomial','FitMethod',fit_method);)
% prv_lme = fitglme(tbl_z,full_prv_mdl,'Distribution','binomial','FitMethod',fit_method)
% % All effects are the same (stim, reward, effort, and stim:reward are
% % significant, but previous reward and effort are not

%% Compute indifference points
% NOPE! This doesn't work... p(Accept) should be 0.5 at that indifference
% point, and it's definitely not.
% ind_on  = log(nthroot(exp(1),par_on(1)));
% pacc_on_ind = (exp(par_on(1)*ind_on) ./ (exp(par_on(1)) + exp(par_on(1)*ind_on)));
% ind_off = log(nthroot(exp(1),par_off(1)));
% pacc_off_ind = (exp(par_off(1)*ind_off) ./ (exp(par_off(1)) + exp(par_off(1)*ind_off)));

%% Average behavior by condition for plotting
% Mean Reward and Effort
efforts = unique(data(:,effort_ix));
stakes  = unique(data(:,stake_ix));
if length(efforts)~=5 || length(stakes)~=5; error('should be 5 conditions!'); end
re_vars = {'Stake','Effort'};
sv_mn             = nan([length(stakes) length(efforts)]);
dec_mn            = nan([length(stakes) length(efforts)]);
dec_mn_on         = nan([length(stakes) length(efforts)]);
dec_mn_off        = nan([length(stakes) length(efforts)]);
pacc_mn          = nan([length(stakes) length(efforts)]);
pacc_mn_on       = nan([length(stakes) length(efforts)]);
pacc_mn_off      = nan([length(stakes) length(efforts)]);
ez_mn             = nan([length(stakes) length(efforts)]);
ez_mn_on          = nan([length(stakes) length(efforts)]);
ez_mn_off         = nan([length(stakes) length(efforts)]);
dec_stake_off_se  = nan([length(stakes) 1]);
dec_stake_on_se   = nan([length(stakes) 1]);
pacc_stake_off_se = nan([length(stakes) 1]);
pacc_stake_on_se  = nan([length(stakes) 1]);
dec_eff_off_se    = nan([length(efforts) 1]);
dec_eff_on_se     = nan([length(efforts) 1]);
pacc_eff_off_se   = nan([length(efforts) 1]);
pacc_eff_on_se    = nan([length(efforts) 1]);
dec_stake_onoff_pval   = nan([length(stakes) 1]);
dec_effort_onoff_pval  = nan([length(stakes) 1]);
pacc_stake_onoff_pval  = nan([length(stakes) 1]);
pacc_effort_onoff_pval = nan([length(stakes) 1]);
for lvl1_i = 1:5
    %-------------- Stake -------------------------
    % OFF stim: Average decision by stake
    stk_idx = data(:,stake_ix)==stakes(lvl1_i);
    cond_idx = stk_idx & off_idx;
    dec_stake_off_se(lvl1_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % OFF stim: Average probability of acceptance by stake
    trl_idx = cond_idx(off_idx);
    pacc_stake_off_se(lvl1_i) = std(pacc_off(trl_idx))./sqrt(sum(trl_idx));
    
    % ON stim: Average decision by stake
    cond_idx = stk_idx & on_idx;
    dec_stake_on_se(lvl1_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % ON stim: Average probability of acceptance by stake
    trl_idx = cond_idx(on_idx);
    pacc_stake_on_se(lvl1_i) = std(pacc_on(trl_idx))./sqrt(sum(trl_idx));
    
    % ON-OFF stim stats: Decision by stake
    [~,dec_stake_onoff_pval(lvl1_i),~,~] = prop_test(...
        [table2array(sum(tbl_z(stk_idx & off_idx,dec_ix))) ...
        table2array(sum(tbl_z(stk_idx & on_idx,dec_ix)))],...
        [sum(stk_idx & off_idx) sum(stk_idx & on_idx)],0);
    pacc_stake_onoff_pval(lvl1_i) = ranksum(pacc_all(stk_idx & off_idx),pacc_all(stk_idx & on_idx));
    % est_dec_lme = fitglme(tbl_z(stk_idx,:),'Decision ~ Effort*Stim',...
    %     'Distribution','binomial','FitMethod','Laplace');
    % [~, est_dec_pval] = fn_GLMM_run_null_LRT_fixed_effects(...
    %     tbl_z(stk_idx,:),est_dec_lme,{'Stim'},'Laplace');
    % dec_stake_onoff_pval(lvl1_i) = est_dec_pval{1}.pValue(2);
    % ON-OFF stim stats: pAccept by stake
    % est_pacc_lme = fitlme(tbl_z(stk_idx,:),'p_accept ~ Effort*Stim','FitMethod','ML');
    % [~, est_pacc_pval] = fn_LMM_run_null_LRT_fixed_effects(...
    %     tbl_z(stk_idx,:),est_pacc_lme,{'Stim'},'ML');
    % pacc_stake_onoff_pval(lvl1_i) = est_pacc_pval{1}.pValue(2);

    %-------------- Effort -------------------------
    % OFF stim: Average decision by effort
    eff_idx = data(:,effort_ix)==efforts(lvl1_i);
    cond_idx = eff_idx & off_idx;
    dec_eff_off_se(lvl1_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % OFF stim: Average probability of acceptance by effort
    trl_idx = cond_idx(off_idx);
    pacc_eff_off_se(lvl1_i) = std(pacc_off(trl_idx))./sqrt(sum(trl_idx));
    
    % ON stim: Average decision by effort
    cond_idx = eff_idx & on_idx;
    dec_eff_on_se(lvl1_i) = std(data(cond_idx,dec_ix))./sqrt(sum(cond_idx));
    % ON stim: Average probability of acceptance by effort
    trl_idx = cond_idx(on_idx);
    pacc_eff_on_se(lvl1_i) = std(pacc_on(trl_idx))./sqrt(sum(trl_idx));

    % ON-OFF stim stats: Decision by effort
    [~,dec_effort_onoff_pval(lvl1_i),~,~] = prop_test(...
        [table2array(sum(tbl_z(eff_idx & off_idx,dec_ix))) ...
        table2array(sum(tbl_z(eff_idx & on_idx,dec_ix)))],...
        [sum(eff_idx & off_idx) sum(eff_idx & on_idx)],0);
    pacc_effort_onoff_pval(lvl1_i) = ranksum(pacc_all(eff_idx & off_idx),pacc_all(eff_idx & on_idx));
    % rst_dec_lme = fitglme(tbl_z(eff_idx,:),'Decision ~ Stake*Stim',...
    %     'Distribution','binomial','FitMethod','Laplace');
    % [~, rst_dec_pval] = fn_GLMM_run_null_LRT_fixed_effects(...
    %     tbl_z(eff_idx,:),rst_dec_lme,{'Stim'},'Laplace');
    % dec_effort_onoff_pval(lvl1_i) = rst_dec_pval{1}.pValue(2);
    % ON-OFF stim stats: pAccept by effort
    % rst_pacc_lme = fitlme(tbl_z(eff_idx,:),'p_accept ~ Stake*Stim','FitMethod','ML');
    % [~, rst_pacc_pval] = fn_LMM_run_null_LRT_fixed_effects(...
    %     tbl_z(eff_idx,:),rst_pacc_lme,{'Stim'},'ML');
    % pacc_effort_onoff_pval(lvl1_i) = rst_pacc_pval{1}.pValue(2);

    for lvl2_i = 1:5
        cond_idx = data(:,effort_ix)==efforts(lvl2_i) & stk_idx;
        % Average Subejctive value
        sv_mn(lvl1_i,lvl2_i) = mean(SV_all(cond_idx));
        % Average over all data
        dec_mn(lvl1_i,lvl2_i) = mean(data(cond_idx,dec_ix));
        pacc_mn(lvl1_i,lvl2_i) = mean(pacc(cond_idx));
        ez_mn(lvl1_i,lvl2_i) = mean(ease_all(cond_idx));
        % Average over OFF data
        trl_idx = cond_idx(off_idx);
        dec_mn_off(lvl1_i,lvl2_i) = mean(data(cond_idx & off_idx,dec_ix));
        pacc_mn_off(lvl1_i,lvl2_i) = mean(pacc_off(trl_idx));
        ez_mn_off(lvl1_i,lvl2_i) = mean(ease_all(cond_idx & off_idx));
        % Average over ON data
        trl_idx = cond_idx(on_idx);
        dec_mn_on(lvl1_i,lvl2_i) = mean(data(cond_idx & on_idx,dec_ix));
        pacc_mn_on(lvl1_i,lvl2_i) = mean(pacc_on(trl_idx));
        ez_mn_on(lvl1_i,lvl2_i) = mean(ease_all(cond_idx & on_idx));
    end
end

% Prob accept ON - OFF
pacc_mn_diff = pacc_mn_on-pacc_mn_off;
dec_mn_diff = dec_mn_on-dec_mn_off;
ez_mn_diff = ez_mn_on-ez_mn_off;

%% Average decision by stim condition
blocks = unique(data(:,blk_ix));
dec_block_off_mn  = nan([length(blocks) 1]);
dec_block_on_mn   = nan([length(blocks) 1]);
for b_ix = 1:length(blocks)
    b_idx = data(:,blk_ix)==blocks(b_ix);
    dec_block_off_mn(b_ix) = mean(data(b_idx & off_idx,dec_ix));
    dec_block_on_mn(b_ix)  = mean(data(b_idx & on_idx,dec_ix));
end
dec_onoff_Xblk_mn = [nanmean(dec_block_off_mn) nanmean(dec_block_on_mn)];
dec_onoff_Xblk_se = [nanstd(dec_block_off_mn)./sqrt(sum(~isnan(dec_block_off_mn))) ...
          nanstd(dec_block_on_mn)./sqrt(sum(~isnan(dec_block_on_mn)))];

% Summary stats for p(Accept)
pacc_onoff_mn = [mean(pacc_off) mean(pacc_on)];
pacc_onoff_se = [nanstd(pacc_off)./sqrt(length(pacc_off)) ...
          nanstd(pacc_on)./sqrt(length(pacc_on))];

%% Bar plot of main effect of stim on work rate (% accept)
fig_name = 'PFC05_stim_allTrials_dec_ONOFF_errbar';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
bars = bar([1 2],diag(dec_onoff_Xblk_mn),'stacked');
bars(1).FaceColor = off_color;
bars(2).FaceColor = on_color;
errorbar([1 2],dec_onoff_Xblk_mn,dec_onoff_Xblk_se,'linewidth',3,'color','k','linestyle','none');
xticklabels(stim_cond);
xticks([1 2]);
ylabel('Work Offeres Accepted (%)');
ylim([0.7 0.81]);
title(['Decision ~ Stimulation (p=' ...
    num2str(dec_comparisons{strcmp(preds,'Stim')}.pValue(2),'%.04f') ')']);
set(gca,'FontSize',18);

if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

fig_name = 'PFC05_stim_allTrials_pacc_ONOFF_errbar';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
bars = bar([1 2],diag(pacc_onoff_mn),'stacked');
bars(1).FaceColor = off_color;
bars(2).FaceColor = on_color;
errorbar([1 2],pacc_onoff_mn,pacc_onoff_se,'linewidth',3,'color','k','linestyle','none');
xticklabels(stim_cond);
xticks([1 2]);
xlabel('PFC Stimulation');
ylabel('Model p(Accept)');
ylim([0.7 0.81]);
title(['Model Probability Accept ~ Stimulation (p=' ...
    num2str(pacc_comparisons{strcmp(preds,'Stim')}.pValue(2),'%.04f') ')']);
set(gca,'FontSize',18);

if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Line plots of Decision and Prob accept by REWARD level
fig_name = 'PFC05_stim_allTrials_dec_stake_ONOFF_line';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(stakes,mean(dec_mn_off,find(strcmp(re_vars,'Effort'))),...
    dec_stake_off_se,'Color',off_color,'linewidth',3);
errorbar(stakes,mean(dec_mn_on,find(strcmp(re_vars,'Effort'))),...
    dec_stake_on_se,'Color',on_color,'linewidth',3);
% for ix = 1:5
%     if dec_stake_onoff_pval(ix)<=0.05
%         y_val = max + se + y_fudge;
%         scatter(stakes(ix),y_val,scat_sz,'k','*');
%     end
% end
xlabel('Reward');
ylabel('Accepted Work Offers (%)');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Decision ~ Stimulation*Reward (p=' ...
    num2str(dec_comparisons{strcmp(preds,'Stake:Stim')}.pValue(2),'%.03f') ')']);
set(gca,'FontSize',18);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/decision_stake_line/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

fig_name = 'PFC05_stim_allTrials_pacc_stake_ONOFF_line';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(stakes,mean(pacc_mn_off,find(strcmp(re_vars,'Effort'))),...
    pacc_stake_off_se,'Color',off_color,'linewidth',3);
errorbar(stakes,mean(pacc_mn_on,find(strcmp(re_vars,'Effort'))),...
    pacc_stake_on_se,'Color',on_color,'linewidth',3);
xlabel('Reward');
ylabel('Model p(Accept)');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Model Probability Accept ~ Stimulation*Reward (p=' ...
    num2str(pacc_comparisons{strcmp(preds,'Stake:Stim')}.pValue(2),'%.04f') ')']);
set(gca,'FontSize',18);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Line plots of Decision and Prob accept by EFFORT level
fig_name = 'PFC05_stim_allTrials_dec_effort_ONOFF_line';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(efforts,mean(dec_mn_off,find(strcmp(re_vars,'Stake'))),...
    dec_eff_off_se,'Color',off_color,'linewidth',3);
errorbar(efforts,mean(dec_mn_on,find(strcmp(re_vars,'Stake'))),...
    dec_eff_on_se,'Color',on_color,'linewidth',3);
xlabel('Effort');
ylabel('Accepted Work Offers (%)');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Decision ~ Stimulation*Effort (p=' ...
    num2str(dec_comparisons{strcmp(preds,'Effort:Stim')}.pValue(2),'%.03f') ')']);
set(gca,'FontSize',18);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/decision_stake_line/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

fig_name = ['PFC05_stim_allTrials_pacc_effort_ONOFF_line'];
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
errorbar(efforts,mean(pacc_mn_off,find(strcmp(re_vars,'Stake'))),...
    pacc_eff_off_se,'Color',off_color,'linewidth',3);
errorbar(efforts,mean(pacc_mn_on,find(strcmp(re_vars,'Stake'))),...
    pacc_eff_on_se,'Color',on_color,'linewidth',3);
xlabel('Effort');
ylabel('Model p(Accept)');
ylim([0 1.1]);
legend({'OFF stimulation','ON stimulation'},'Location','best');
title(['Model Probability Accept ~ Stimulation*Effort (p=' ...
    num2str(pacc_comparisons{strcmp(preds,'Effort:Stim')}.pValue(2),'%.04f') ')']);
set(gca,'FontSize',18);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Line plot for SV vs. decision function
fig_name = 'PFC05_stim_dec_fn_ONOFF_line';
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
off_line = line(SV_off(SV_off_sort_idx),pacc_off(SV_off_sort_idx),'Color',off_color,'LineWidth',3);
% scatter(SV_off(SV_off_sort_idx),pacc_off(SV_off_sort_idx),scat_sz,off_color);
on_line = line(SV_on(SV_on_sort_idx),pacc_on(SV_on_sort_idx),'Color',on_color,'LineWidth',3);
% scatter(SV_on(SV_on_sort_idx),pacc_on(SV_on_sort_idx),scat_sz,on_color);
%     line(mdl_x,mdl_y,'Color','k');
line([-5 15],[0.5 0.5],'Color','k','LineStyle','--');
xlim([-5 15]);
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
end

sv_vals  = -4.5:0.01:13;
pacc_vals_on  = (exp(par_on(1)*sv_vals)) ./ (exp(par_on(1)) + exp(par_on(1)*sv_vals));
pacc_vals_off = (exp(par_off(1)*sv_vals)) ./ (exp(par_off(1)) + exp(par_off(1)*sv_vals));

% fig_name = 'PFC05_stim_dec_fn_ONOFF_line_fine';
% figure('Name',fig_name,'units','norm','outerposition',[0 0 0.3 0.5]); hold on;
% off_line = line(sv_vals,pacc_vals_off,'Color',off_color,'LineWidth',3);
% on_line = line(sv_vals,pacc_vals_on,'Color',on_color,'LineWidth',3);
% line([-5 15],[0.5 0.5],'Color','k','LineStyle','--');
% xlim([-5 15]);
% xlabel('Subjective Value'); ylabel('Model p(Accept)');
% title('Subjective Value Decision Function');
% legend([off_line,on_line],{'OFF stimulation','ON stimulation'},'Location','best');
% set(gca,'FontSize',font_sz);

%% 3D bar plot of subjective value landscape
fig_name = 'PFC05_stim_allTrials_sv_bar3';
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,sv_mn,font_sz);
zlabel('Subjective Value');
title('PFC05 Subjective Value Landscape: All trials');

% Shade bar according to z-axis
for b = 1:length(bars)
    zdata = bars(b).ZData;
    bars(b).CData = zdata;
    bars(b).FaceColor = 'interp';
end
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/sv_bar3/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% 3-D plot of decision as function of reward and effort
% Total behavior
fig_name = ['PFC05_stim_allTrials_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,dec_mn,font_sz);
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
[bars] = fn_3d_decision_landscape(stakes,efforts,dec_mn_diff,font_sz);
zlabel('% Accept');
title('PFC05 Decision Landscape: ON-OFF');
rb_cmap = redblue();
colormap(rb_cmap);
caxis([-max(abs(dec_mn_diff(:))) max(abs(dec_mn_diff(:)))]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% OFF stim
fig_name = ['PFC05_stim_OFF_p_dec_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,dec_mn_off,font_sz);
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
[bars] = fn_3d_decision_landscape(stakes,efforts,dec_mn_on,font_sz);
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
fig_name = ['PFC05_stim_allTrials_pacc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,pacc_mn,font_sz);

title('PFC05 Decision Model Landscape: All trials');
caxis([0 1]);
if save_fig
    fig_dir   = [prj_dir 'results/bhv/PFC05_stim/paccept_bar3/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON-OFF behavior
fig_name = ['PFC05_stim_ON-OFF_pacc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,pacc_mn_diff,font_sz);
title('PFC05 Decision Model Landscape: ON-OFF');
rb_cmap = redblue();
colormap(rb_cmap);
caxis([-max(abs(pacc_mn_diff(:))) max(abs(pacc_mn_diff(:)))]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% 2D matrix of ON-OFF behavior
fig_name = ['PFC05_stim_ON-OFF_pacc_mat'];
figure('Name',fig_name);
imagesc(1:length(stakes),1:length(efforts),pacc_mn_diff);
set(gca,'YDir','normal');
yticks(1:length(efforts));
yticklabels(stakes);
xticks(1:length(stakes));
xlabel('Effort (% Max)');
xticklabels(efforts*100);
cbar = colorbar;
ylabel(cbar,'Diff. Prob. Accept (ON-OFF)','Rotation',270);
rb_cmap = redblue();
colormap(rb_cmap);
caxis([-max(abs(pacc_mn_diff(:))) max(abs(pacc_mn_diff(:)))]);
set(gca,'FontSize',font_sz);
title('PFC05 Decision Model Landscape: ON-OFF');
cbar.Label.Position(1) = 4;
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% OFF stim
fig_name = ['PFC05_stim_OFF_pacc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,pacc_mn_off,font_sz);
title('PFC05 Decision Model Landscape: OFF Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

% ON stim
fig_name = ['PFC05_stim_ON_pacc_bar3'];
figure('Name',fig_name);
[bars] = fn_3d_decision_landscape(stakes,efforts,pacc_mn_on,font_sz);
title('PFC05 Decision Model Landscape: ON Stimulation');
caxis([0 1]);
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Model behavior in difficult decision range
low_ez_idx = ease_all<median(ease_all);
% Redo zscore with only relevant data
zReg_lowez = [zscore(data(low_ez_idx,effort_ix)) zscore(data(low_ez_idx,stake_ix)) data(low_ez_idx,3:7) zscore(paccept_all(low_ez_idx))];
zT_lowez = array2table(zReg_midrew, 'VariableNames', col_names);

full_mdl = 'p_accept ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block) + (1|TrialCum) + (1|BlkCum)';
zmdl_pacc_lowez = fitlme(zT_lowez,full_mdl)

% Compare to a model with only ease
zReg_lowez = [zscore(data(low_ez_idx,effort_ix)) zscore(data(low_ez_idx,stake_ix)) data(low_ez_idx,3:7) zscore(paccept_all(low_ez_idx))];
zT_lowez = array2table(zReg_midrew, 'VariableNames', col_names);

full_mdl = 'p_accept ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block) + (1|TrialCum) + (1|BlkCum)';
zmdl_pacc_lowez = fitlme(zT_lowez,full_mdl)

% Plot correlation between p accept ON-OFF and decision difficuty

figure;
scatter(sv_mn(:),pacc_mn_diff(:));
xlabel('SV');
ylabel('Prob. Accept ON-OFF');

%% Model behavior in middle reward range
midstake_idx = data(:,stake_ix)==4 | data(:,stake_ix)==7;
% Redo zscore with only relevant data
zReg_midrew = [zscore(data(midstake_idx,effort_ix)) zscore(data(midstake_idx,stake_ix)) data(midstake_idx,3:7) zscore(paccept_all(midstake_idx))];
zReg_midrew_ease = [zscore(data(midstake_idx,effort_ix)) zscore(data(midstake_idx,stake_ix)) data(midstake_idx,3:7) zscore(paccept_all(midstake_idx)) zscore(ease_all(midstake_idx))];
% T = array2table(data, 'VariableNames', col_names);
zT_midrew = array2table(zReg_midrew, 'VariableNames', col_names);
zT_midrew_ease = array2table(zReg_midrew_ease, 'VariableNames', [col_names 'ease']);

full_mdl = 'p_accept ~ Effort*Stake*Stim - Effort:Stake:Stim + (1|Block) + (1|TrialCum) + (1|BlkCum)';
zmdl_pacc_stakeeff47 = fitlme(zT_midrew,full_mdl)
modelspec_ez = 'p_accept ~ ease*Stim + (1|Block) + (1|TrialCum) + (1|BlkCum)';
zmdl_pacc_ez47 = fitlme(zT_midrew_ease,modelspec_ez)

% Check mid-range stakes where difference is most pronounced
% stake_val = 4;
% stake_idx = data(:,stake_ix)==stake_val;
% modelspec = 'p_accept ~ Stim*Effort + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
% zmdl_pacc_stake4 = fitlme(zT(stake_idx,:),modelspec)
% 
% stake_val = 7;
% stake_idx = data(:,stake_ix)==stake_val;
% modelspec = 'p_accept ~ Stim*Effort + (1|Block)';% + (1|TrialCum) + (1|BlkCum)';
% zmdl_pacc_stake7 = fitlme(zT(stake_idx,:),modelspec)

%% define function for 3D bar plot
function [bars] = fn_3d_decision_landscape(stakes,efforts,zvals,font_sz)

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
%     idxrep = repmat(1:size(pacc_mn_stim,1),6,1);  %magic 6
%     for col_ix = 1:size(pacc_mn_stim,2)
%         Z = repmat(pacc_mn_stim(idxrep(:),col_ix,p),1,4);   %magic 4
%         set(bars(col_ix), 'CData',Z);
%     end
end

