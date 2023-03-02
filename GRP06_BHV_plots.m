%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

save_fig  = 1;
fig_ftype = 'svg';%'png';

addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

%% Load data
bhvs       = cell(size(SBJs));
mdls       = cell(size(SBJs));
for s = 1:length(SBJs)
    % Load behavior
    load([prj_dir 'data/' SBJs{s} '/' SBJs{s} '_stim_preproc.mat'],'sbj_data');
    bhvs{s} = sbj_data.bhv;
    mdls{s} = sbj_data.mdl;
end

%% Behavioral model fitting
%   SV_physical(t) = R(t)?k*E(t)^2
%   softmax: exp(B * Q1) / exp(B * Q1) + exp(B * Q2)
%   here, Q1 is SV and Q2 is nothing (reject choice, no alternative, so just B term)
% decisionfun=@(p) norm( (exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))) ./ ...
%     (exp(p(1)) + exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))))) - bhv.decision);
% [par, fit]=fminsearch(decisionfun, [1,1]);
% 
% SV_fn    = @(k) bhv.stake-(k*(bhv.effort).^2);
% EFF_fn   = @(k) (k*(bhv.effort).^2);
% bhv.SV   = SV_fn(par(2));
% bhv.EFFs = EFF_fn(par(2));
% bhv.p_accept = (exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2))) ./...
%     (exp(par(1)) + exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2)))));

%% Create variables
% Mean Reward and Effort
efforts = unique(bhvs{1}.effort);
stakes  = unique(bhvs{1}.stake);
p_acc_mn = nan([length(stakes) length(efforts) length(SBJs)]);
for s = 1:length(SBJs)
    if length(unique(bhvs{s}.effort))~=5 || length(unique(bhvs{s}.stake))~=5
        error([SBJs{s} ' is missing conditions for stake or effort!']);
    end
    bhvs{s}.SV_stake_mn  = nan([5 1]);
    bhvs{s}.SV_stake_se  = nan([5 1]);
    bhvs{s}.SV_effort_mn = nan([5 1]);
    bhvs{s}.SV_effort_se = nan([5 1]);
    for stk_i = 1:5
        bhvs{s}.SV_effort_mn(stk_i) = mean(bhvs{s}.SV(bhvs{s}.effort==efforts(stk_i)));
        bhvs{s}.SV_effort_se(stk_i) = std(bhvs{s}.SV(bhvs{s}.effort==efforts(stk_i)))./sqrt(sum(bhvs{s}.effort==efforts(stk_i)));
        bhvs{s}.SV_stake_mn(stk_i) = mean(bhvs{s}.SV(bhvs{s}.stake==stakes(stk_i)));
        bhvs{s}.SV_stake_se(stk_i) = std(bhvs{s}.SV(bhvs{s}.stake==stakes(stk_i)))./sqrt(sum(bhvs{s}.stake==stakes(stk_i)));
        for eff_i = 1:5
            p_acc_mn(stk_i,eff_i,s) = mean(bhvs{s}.p_accept(bhvs{s}.effort==efforts(eff_i) & bhvs{s}.stake==stakes(stk_i)));
        end
    end
end

p_acc_mn_grp = mean(p_acc_mn,3);

%% 3-D plot of acceptance as function of reward and effort
% Subject level
for s = 1:length(SBJs)
    fig_name = [SBJs{s} '_p_acc_bar3'];
    figure('Name',fig_name);
    bars = bar3(squeeze(p_acc_mn(:,:,s)));
    ylabel('Reward');
    yticklabels(stakes);
    xlabel('Effort (% Max Force)');
    xticklabels(efforts*100);
    zlabel('% Accept');
    title([SBJs{s} ' Decision Landscape']);
    set(gca,'FontSize',16);
    view([135 30]);
    
    % Shade bar according to z-axis
    for b = 1:length(bars)
        zdata = bars(b).ZData;
        bars(b).CData = zdata;
        bars(b).FaceColor = 'interp';
    end
%     idxrep = repmat(1:size(p_acc_mn,1),6,1);  %magic 6
%     for col_ix = 1:size(p_acc_mn,2)
%         Z = repmat(p_acc_mn(idxrep(:),col_ix),1,4);   %magic 4
%         set(bars(col_ix), 'CData',Z);
%     end
    
    if save_fig
        fig_dir   = [prj_dir 'results/bhv/p_accept_bar3/'];
        if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

% Group level
fig_name = ['GRP_p_acc_bar3'];
figure('Name',fig_name);
bars = bar3(p_acc_mn_grp);
ylabel('Reward');
yticklabels(stakes);
xlabel('Effort (% Max)');
xticklabels(efforts*100);
zlabel('% Accept');
title('Group Decision Landscape');
set(gca,'FontSize',16);
view([135 30]);

% Shade bar according to z-axis
for b = 1:length(bars)
    zdata = bars(b).ZData;
    bars(b).CData = zdata;
    bars(b).FaceColor = 'interp';
end
% idxrep = repmat(1:size(p_acc_mn_grp,1),6,1);  %magic 6
% for col_ix = 1:size(p_acc_mn_grp,2)
%   Z = repmat(p_acc_mn_grp(idxrep(:),col_ix),1,4);   %magic 4
%   set(bars(col_ix), 'CData',Z);
% end

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot correlation between behavioral predictors
scat_sz = 40;

pred = {'rt','effort','stake','decision','SV','EFFs','p_accept','absSV','dec_diff'};
pred_prv = strcat(pred,'_prv');
pred_all = [pred pred_prv];
prv_idx = contains(pred_all,'_prv');
pred_mat = cell(size(SBJs));
for s = 1:length(SBJs)
    % Build matrix
    pred_mat{s} = nan([length(bhvs{s}.rt) length(pred_all)]);
    for p = 1:length(pred_all)
        pred_mat{s}(:,p) = bhvs{s}.(pred_all{p});
    end
    
    % Compute correlations
    [pred_all_corr, pred_all_pval] = corrcoef(pred_mat{s},'Rows','complete');
    pred_all_sig = pred_all_pval<=0.05;
    
    % Plot current correlations
    fig_name = [SBJs{s} '_bhv_pred_corr_mat'];
    figure('Name',fig_name); hold on;
    imagesc(pred_all_corr(prv_idx==0,prv_idx==0));
    set(gca,'XTick',1:length(pred));
    set(gca,'XTickLabels',pred);
    set(gca,'XTickLabelRotation',45);
    set(gca,'XLim',[0.5 length(pred)+0.5]);
    set(gca,'YTick',1:length(pred));
    set(gca,'YTickLabels',pred);
    set(gca,'YTickLabelRotation',45);
    set(gca,'YLim',[0.5 length(pred)+0.5]);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',16);
    colorbar;
    caxis([-1 1]);
    title('Behavioral Predictor Correlations');
    for p1 = 1:length(pred)
        for p2 = 1:length(pred)
            if pred_all_sig(p1,p2); scatter(p1,p2,scat_sz,'k','*'); end
        end
    end
    if save_fig
        fig_dir   = [prj_dir 'results/bhv/pred_corr/'];
        if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
    
    % Plot current vs. previous correlations
    fig_name = [SBJs{s} '_bhv_pred_curVsPrv_corr_mat'];
    figure('Name',fig_name); hold on;
    imagesc(pred_all_corr(prv_idx==0,prv_idx==1));
    set(gca,'XTick',1:length(pred));
    set(gca,'XTickLabels',pred);
    set(gca,'XTickLabelRotation',45);
    set(gca,'XLim',[0.5 length(pred)+0.5]);
    set(gca,'YTick',1:length(pred));
    set(gca,'YTickLabels',pred_prv);
    set(gca,'YTickLabelRotation',45);
    set(gca,'YLim',[0.5 length(pred)+0.5]);
    set(gca,'TickLabelInterpreter','none');
    set(gca,'FontSize',16);
    colorbar;
    title('Behavioral Predictor Correlations: Current vs. Previous');
%     caxis([-1 1]);
    for p1 = 1:length(pred)
        for p2 = 1+length(pred):length(pred_all)
            if pred_all_sig(p1,p2); scatter(p1,p2-length(pred_all),scat_sz,'k','*'); end
        end
    end
    
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot effort vs. EFF
fig_name = 'GRP_effort_vs_EFFs';
figure('Name',fig_name);%,'units','norm','OuterPosition',[0 0 0.5 0.6]);
hold on;
for s = 1:length(SBJs)
    [~,eff_sort_idx] = sort(bhvs{s}.effort);
    line(bhvs{s}.effort(eff_sort_idx),bhvs{s}.EFFs(eff_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',2);
    xlabel('Effort');
    ylabel('Subjective Effort');
    title('Subjective Effort Function');
    set(gca,'FontSize',16);
end
legend(SBJs,'Location','best');

if save_fig
    fig_dir   = [prj_dir 'results/bhv/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Characterize behavior: Reward, subejctive ffort, SV decision fn, abs(SV)
scat_sz = 40;
mdl_x = -15:0.001:15;
fig_name = ['GRP_bhv_model_scatter'];
figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.5 0.6]);
for s = 1:length(SBJs)
%     mdl_y = mdls{s}.par(1) + (mdl_x * mdls{s}.par(2));
%     mdl_y = 1 ./ (1+exp(-mdl_y));
%     mdl_y_vals = (exp(mdls{s}.par(1)*mdl_x) ./ ...
%         (exp(mdls{s}.par(1)) + exp(mdls{s}.par(1)*mdl_x)));
%     betas = glmfit(bhvs{s}.SV,bhvs{s}.p_accept,'binomial','link','logit');
%     sigparam = sigm_fit(bhvs{s}.SV,bhvs{s}.p_accept);
    [~,SV_sort_idx] = sort(bhvs{s}.SV);
    
    subplot(2,2,1); hold on;
%     scatter(bhvs{s}.stake,bhvs{s}.SV,scat_sz,sbj_colors(s,:));
    errorbar(stakes,bhvs{s}.SV_stake_mn,bhvs{s}.SV_stake_se,...
        'Color',sbj_colors(s,:),'LineWidth',2);
    xlabel('Reward'); ylabel('SV');
    title('Reward vs. Subjective Value');
    set(gca,'FontSize',16);
    legend(SBJs,'Location','northwest');
    
    subplot(2,2,2); hold on;
%     scatter(bhvs{s}.effort,bhvs{s}.SV,scat_sz,sbj_colors(s,:));
    errorbar(efforts,bhvs{s}.SV_effort_mn,bhvs{s}.SV_effort_se,...
        'Color',sbj_colors(s,:),'LineWidth',2);
    xlabel('Effort - Proportion max'); ylabel('SV');
    title('Objective Effort vs. Subjective Value');
    set(gca,'FontSize',16);
    
    subplot(2,2,3); hold on;
    line(bhvs{s}.SV(SV_sort_idx),bhvs{s}.p_accept(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',2);
    scatter(bhvs{s}.SV,bhvs{s}.p_accept,scat_sz,sbj_colors(s,:));
%     line(mdl_x,mdl_y,'Color','k');
    line([0 0],[0 1],'Color','k','LineStyle','--');
    xlabel('SV'); ylabel('Probabilty of acceptance');
    title('Subjective Value Decision Function');
    set(gca,'FontSize',16);

    subplot(2,2,4); hold on;
    line(abs(bhvs{s}.SV(SV_sort_idx)-median(bhvs{s}.SV)),bhvs{s}.p_accept(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',2);
    scatter(abs(bhvs{s}.SV-median(bhvs{s}.SV)),bhvs{s}.p_accept,scat_sz,sbj_colors(s,:));
%     line(mdl_x,mdl_y,'Color','k');
    xlabel('abs(SV)'); ylabel('Probabilty of acceptance');
    title('Absolute Subjective Value');
    set(gca,'FontSize',16);
end

if save_fig
    fig_dir   = [prj_dir 'results/bhv/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Clean plot of subjective value decision function
fig_name = ['GRP_bhv_model_SV_decision_fn'];
figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.3 0.4]);
hold on;
for s = 1:length(SBJs)
    [~,SV_sort_idx] = sort(bhvs{s}.SV);
    
    line(bhvs{s}.SV(SV_sort_idx),bhvs{s}.p_accept(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',3);
%     scatter(bhvs{s}.SV,bhvs{s}.p_accept,scat_sz,sbj_colors(s,:));
%     line(mdl_x,mdl_y,'Color','k');
    line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
    xlabel('Subjective Value'); ylabel('Probabilty Accept');
    title('Subjective Value Decision Function');
    set(gca,'FontSize',18);
end

if save_fig
    fig_dir   = [prj_dir 'results/bhv/'];
    if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Split decision by input variables
for s=1:4
    fig_name = [SBJs{s} '_bhv_model_hist_split_dec'];
    figure('Name',fig_name);
    subplot(2,2,1); hold on;
    histogram(bhvs{s}.SV(bhvs{s}.decision==0),'FaceColor','r','FaceAlpha',0.3);
    histogram(bhvs{s}.SV(bhvs{s}.decision==1),'FaceColor','b','FaceAlpha',0.3);
    xlabel('Current SV');
    title(SBJs{s});
    legend('Reject','Accept');
    set(gca,'FontSize',16);
    
    subplot(2,2,2); hold on;
    histogram(bhvs{s}.SV_prv(bhvs{s}.decision==0),'FaceColor','r','FaceAlpha',0.3);
    histogram(bhvs{s}.SV_prv(bhvs{s}.decision==1),'FaceColor','b','FaceAlpha',0.3);
    xlabel('Previous Trial SV');
    legend('Reject','Accept');
    set(gca,'FontSize',16);
    
    subplot(2,2,3); hold on;
    histogram(bhvs{s}.stake(bhvs{s}.decision==0),'FaceColor','r','FaceAlpha',0.3);
    histogram(bhvs{s}.stake(bhvs{s}.decision==1),'FaceColor','b','FaceAlpha',0.3);
    xlabel('Stake/reward');
    legend('Reject','Accept');
    set(gca,'FontSize',16);
    
    subplot(2,2,4); hold on;
    histogram(bhvs{s}.EFFs(bhvs{s}.decision==0),'FaceColor','r','FaceAlpha',0.3);
    histogram(bhvs{s}.EFFs(bhvs{s}.decision==1),'FaceColor','b','FaceAlpha',0.3);
    xlabel('EFFs');
    legend('Reject','Accept');
    set(gca,'FontSize',16);
    
    if save_fig
        fig_dir   = [prj_dir 'results/bhv/'];
        if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Plot a regression line %
% SVtmp = bhvs{s}.SV;
% 
% figure; scatter(SVtmp,ProbAccept,'.')
% 
% sigparam=sigm_fit(SVtmp,ProbAccept,1)
% 
% [p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
% PotValues=[-10:0.1:10];
% Pout=polyval(p,PotValues);
% 
% figure; plot(PotValues,Pout)
% 
% % Test behavioral modelling.
% figure;
% subplot(2,1,1);
% scatter(bhvs{s}.stake, bhvs{s}.SV);
% title('Reward')
% subplot(2,1,2);
% scatter(bhvs{s}.effort, bhvs{s}.SV);
% title('Effort');
% 
