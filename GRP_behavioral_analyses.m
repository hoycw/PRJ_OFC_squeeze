%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
SBJs         = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
sbj_bg_roi   = {'GPi','STN','GPi','STN'};

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);

%% Load data
bhvs       = cell(size(SBJs));
mdls       = cell(size(SBJs));
for s = 1:length(SBJs)
    % Load behavior
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    proc_fname = [sbj_dir SBJs{s} '_stim_preproc.mat'];
    fprintf('Loading %s\n',proc_fname);
    load(proc_fname,'tfr');
    load([sbj_dir SBJs{s} '_stim_preproc.mat']);
    bhvs{s} = sbj_data.bhv;
    mdls{s} = sbj_data.mdl;
end

%% Model fitting
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

%% Are decisions related to SV?
mdl_x = -15:0.001:15;
for s = 1:length(SBJs)
%     mdl_y = mdls{s}.par(1) + (mdl_x * mdls{s}.par(2));
%     mdl_y = 1 ./ (1+exp(-mdl_y));
%     mdl_y_vals = (exp(mdls{s}.par(1)*mdl_x) ./ ...
%         (exp(mdls{s}.par(1)) + exp(mdls{s}.par(1)*mdl_x)));
    betas = glmfit(bhvs{s}.SV,bhvs{s}.p_accept,'binomial','link','logit');
    sigparam = sigm_fit(bhvs{s}.SV,bhvs{s}.p_accept);
    
    fig_name = [SBJs{s} '_bhv_model_scatter'];
    figure('Name',fig_name);
    subplot(3,1,1);
    scatter(bhvs{s}.stake,bhvs{s}.SV);
    xlabel('Reward'); ylabel('SV');
    title('Reward vs. Subjective Value');
    set(gca,'FontSize',16);
    subplot(3,1,2);
    scatter(bhvs{s}.effort,bhvs{s}.SV);
    xlabel('Effort - Proportion max'); ylabel('SV');
    title('Objective Effort vs. Subjective Value');
    set(gca,'FontSize',16);
    subplot(3,1,3);
    scatter(bhvs{s}.SV,bhvs{s}.p_accept);
    line(mdl_x,mdl_y,'Color','k');
    xlabel('SV'); ylabel('Probabilty of acceptance');
    title('Subjective Value Decision Function');
    set(gca,'FontSize',16);
end

%%
figure; hold on;
histogram(bhvs{s}.SV(bhvs{s}.decision==0),'FaceColor','r','FaceAlpha',0.3);
histogram(bhvs{s}.SV(bhvs{s}.decision==1),'FaceColor','b','FaceAlpha',0.3);
legend('Reject','Accept');

%% Plot a regression line %
SVtmp = bhvs{s}.SV;

figure; scatter(SVtmp,ProbAccept,'.')

sigparam=sigm_fit(SVtmp,ProbAccept,1)

[p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
PotValues=[-10:0.1:10];
Pout=polyval(p,PotValues);

figure; plot(PotValues,Pout)

% Test behavioral modelling.
figure;
subplot(2,1,1);
scatter(bhvs{s}.stake, bhvs{s}.SV);
title('Reward')
subplot(2,1,2);
scatter(bhvs{s}.effort, bhvs{s}.SV);
title('Effort');

