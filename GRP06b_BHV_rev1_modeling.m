%% Mixed modelling on the Processed data %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/sigm_fit/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults
close all
clear all

%%
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'; stat_id = 'S5t15_bhv0_nrl0_out3';

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
    load([prj_dir 'data/' SBJs{s} '/' SBJs{s} '_stim_preproc_osr.mat'],'sbj_data');
    bhvs{s} = sbj_data.bhv;
    mdls{s} = sbj_data.mdl;
end

%% Load group model tables
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
tbl = readtable(table_all_fname);

%% Behavioral model fitting
% Original model:
%   SV_physical(t) = R(t) - k*E(t)^2
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

%% Revision 1 model to integrate effect of previous reward:
%   1) Add personal reward term:
%       SV = k1*R - k2*E^2
%   2) Add previous reward interaction:
%       SV = k1*(R - R_prv) - k2*E^2

%% Original model
options = optimset('PlotFcns',@optimplotfval);%,'MaxFunEvals',5000,'MaxIter',5000);
orig_par = nan(length(SBJs),2);
orig_fit = nan(size(SBJs));
orig_flg = nan(size(SBJs));
for s = 1:length(SBJs)
    s_idx = tbl.sbj_n==s;
    dec_fn=@(p) norm( (exp(p(1)*(tbl.reward_cur(s_idx)-(p(2)*(tbl.effort_cur(s_idx)).^2))) ./ ...
        (exp(p(1)) + exp(p(1)*(tbl.reward_cur(s_idx)-(p(2)*(tbl.effort_cur(s_idx)).^2))))) - tbl.decision_cur(s_idx));
    [orig_par(s,:), orig_fit(s), orig_flg(s)]=fminsearch(dec_fn, [1,1], options);
    orig_SV_fn    = @(k) tbl.reward_cur(s_idx)-(k*(tbl.effort_cur(s_idx)).^2);
    orig_EFF_fn   = @(k) (k*(tbl.effort_cur(s_idx)).^2);
    SV{s}.orig   = orig_SV_fn(orig_par(s,2));
    EFFs{s}.orig = orig_EFF_fn(orig_par(s,2));
    pacc{s}.orig = (exp(orig_par(s,1)*(tbl.reward_cur(s_idx)-(orig_par(s,2)*(tbl.effort_cur(s_idx)).^2))) ./...
        (exp(orig_par(s,1)) + exp(orig_par(s,1)*(tbl.reward_cur(s_idx)-(orig_par(s,2)*(tbl.effort_cur(s_idx)).^2)))));
    % pause;
    tbl.pacc_orig(s_idx) = pacc{s}.orig;
end
% PFC05 doesn't converge to reasonable values, so maybe this isn't a good
% move...

fig_name = 'GRP_bhv_model_orig_SV_decision_fn';
figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.3 0.4]);
hold on;
sbj_idx = [1 2 3 4];
for s = sbj_idx
    [~,SV_sort_idx] = sort(SV{s}.orig);
    
    sv_lines(s) = line(SV{s}.orig(SV_sort_idx),pacc{s}.orig(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',3);
%     scatter(bhvs{s}.SV,bhvs{s}.p_accept,scat_sz,sbj_colors(s,:));
%     line(mdl_x,mdl_y,'Color','k');
end
line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
xlabel('Subjective Value'); ylabel('Probabilty Accept');
title('Original SV Decision Function');
legend(sv_lines,SBJs(sbj_idx));
set(gca,'FontSize',18);

%% Adding individualized reward term
options = optimset('PlotFcns',@optimplotfval);%,'MaxFunEvals',5000,'MaxIter',5000);
irew_par = nan(length(SBJs),3);
irew_fit = nan(size(SBJs));
irew_flg = nan(size(SBJs));
for s = 1:length(SBJs)
    s_idx = tbl.sbj_n==s;
    irew_dec_fn=@(p) norm( (exp(p(1)*((p(2)*tbl.reward_cur(s_idx))-(p(3)*(tbl.effort_cur(s_idx)).^2))) ./ ...
        (exp(p(1)) + exp(p(1)*((p(2)*tbl.reward_cur(s_idx))-(p(3)*(tbl.effort_cur(s_idx)).^2))))) - tbl.decision_cur(s_idx));
    [irew_par(s,:), irew_fit(s), irew_flg(s)]=fminsearch(irew_dec_fn, [1,1,1], options);
    irew_SV_fn    = @(k) (k(1)*tbl.reward_cur(s_idx))-(k(2)*(tbl.effort_cur(s_idx)).^2);
    irew_EFF_fn   = @(k) (k*(tbl.effort_cur(s_idx)).^2);
    SV{s}.irew   = irew_SV_fn(irew_par(s,2:3));
    EFFs{s}.irew = irew_EFF_fn(irew_par(s,3));
    pacc{s}.irew = (exp(irew_par(s,1)*((irew_par(s,3)*tbl.reward_cur(s_idx))-(irew_par(s,2)*(tbl.effort_cur(s_idx)).^2))) ./...
        (exp(irew_par(s,1)) + exp(irew_par(s,1)*(tbl.reward_cur(s_idx)-(irew_par(s,2)*(tbl.effort_cur(s_idx)).^2)))));
    % pause;
end
% PFC05 doesn't converge to reasonable values, so maybe this isn't a good
% move...

fig_name = 'GRP_bhv_model_indRew_SV_decision_fn';
figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.3 0.4]);
hold on;
sbj_idx = [1 2 3 4];
for s = sbj_idx
    [~,SV_sort_idx] = sort(SV{s}.irew);
    
    sv_line(s) = line(SV{s}.irew(SV_sort_idx),pacc{s}.irew(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',3);
%     scatter(bhvs{s}.SV,bhvs{s}.p_accept,scat_sz,sbj_colors(s,:));
%     line(mdl_x,mdl_y,'Color','k');
end
line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
xlabel('Subjective Value'); ylabel('Probabilty Accept');
title('Individualized Reward SV Decision Function');
legend(sv_lines,SBJs(sbj_idx));
set(gca,'FontSize',18);

% if save_fig
%     fig_dir   = [prj_dir 'results/bhv/'];
%     if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end

%% Add previous reward interactions
options = optimset('PlotFcns',@optimplotfval);%,'MaxFunEvals',5000,'MaxIter',5000);
sprew_par = nan(length(SBJs),2);
sprew_fit = nan(size(SBJs));
sprew_flg = nan(size(SBJs));
pacc_sprew = nan(size(tbl.pAccept_cur));
for s = 1:length(SBJs)
    s_idx = tbl.sbj_n==s & ~isnan(tbl.reward_prv);
    sprew_dec_fn=@(p) norm( (exp(p(1)*((tbl.reward_cur(s_idx)-tbl.reward_prv(s_idx))-...
        (p(2)*(tbl.effort_cur(s_idx)).^2))) ./ ...
        (exp(p(1)) + exp(p(1)*((tbl.reward_cur(s_idx)-tbl.reward_prv(s_idx))-...
        (p(2)*(tbl.effort_cur(s_idx)).^2))))) - tbl.decision_cur(s_idx));
    [sprew_par(s,:), sprew_fit(s), sprew_flg(s)]=fminsearch(sprew_dec_fn, [1,1], options);
    sprew_SV_fn    = @(k) tbl.reward_cur(s_idx)-(k*(tbl.effort_cur(s_idx)).^2);
    sprew_EFF_fn   = @(k) (k*(tbl.effort_cur(s_idx)).^2);
    SV{s}.sprew   = sprew_SV_fn(sprew_par(s,2));
    EFFs{s}.sprew = sprew_EFF_fn(sprew_par(s,2));
    pacc{s}.sprew = (exp(sprew_par(s,1)*(tbl.reward_cur(s_idx)-(sprew_par(s,2)*(tbl.effort_cur(s_idx)).^2))) ./...
        (exp(sprew_par(s,1)) + exp(sprew_par(s,1)*(tbl.reward_cur(s_idx)-(sprew_par(s,2)*(tbl.effort_cur(s_idx)).^2)))));
    % pause;
    pacc_sprew(s_idx) = pacc{s}.sprew;
end
tbl.pacc_sprew = pacc_sprew;

%%
fig_name = 'GRP_bhv_model_subPrvRew_SV_decision_fn';
figure('Name',fig_name,'units','norm','OuterPosition',[0 0 0.5 0.4]);
subplot(1,3,1); hold on;
sbj_idx = [1 2 3 4];
subplot(1,3,1); hold on;
for s = sbj_idx
    [~,SV_sort_idx] = sort(SV{s}.sprew);
    sv_lines(s) = plot(SV{s}.sprew(SV_sort_idx),pacc{s}.sprew(SV_sort_idx),...
        'Color',sbj_colors(s,:),'LineWidth',3);
    scatter(SV{s}.sprew,pacc{s}.sprew,75,sbj_colors(s,:));
end
line(xlim,[0.5 0.5],'Color','k','LineStyle','--');
xlabel('Subjective Value'); ylabel('Probabilty Accept');
title('Subtract Prev Reward SV Decision Function');
legend(sv_lines,SBJs(sbj_idx));
set(gca,'FontSize',18);

subplot(1,3,2); hold on;
fn_plot_LMM_quantile_line_interaction(tbl(~isnan(tbl.reward_prv),:),'reward_prv','reward_cur','pAccept_cur',2,5);
subplot(1,3,3); hold on;
fn_plot_LMM_quantile_line_interaction(tbl(~isnan(tbl.reward_prv),:),'reward_prv','reward_cur','pacc_sprew',2,5);

figure;
for s = 1:length(SBJs)
    subplot(2,2,s); hold on;
    fn_plot_LMM_quantile_line_interaction(tbl(~isnan(tbl.reward_prv) & tbl.sbj_n==s,:),'reward_prv','reward_cur','pAccept_cur',2,5);
    title(SBJs{s});
end
figure;
for s = 1:length(SBJs)
    subplot(2,2,s); hold on;
    fn_plot_LMM_quantile_line_interaction(tbl(~isnan(tbl.reward_prv) & tbl.sbj_n==s,:),'reward_prv','reward_cur','pacc_sprew',2,5);
    title(SBJs{s});
end

% if save_fig
%     fig_dir   = [prj_dir 'results/bhv/'];
%     if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
%     fig_fname = [fig_dir fig_name '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end
