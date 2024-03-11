%% Plot time-resolved LMM coefficients
clear all
close all

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Analysis parameters:
an_id = 'TFRmth_S1t2_madS8t0_f2t40_osr'; stat_id = 'S0t2ts_wl2s025_bhvz_nrl0_out3main';
pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
predictors  = {'reward_cur','effortS_cur','reward_prv','effortS_prv'};
% for p = 1:length(predictors)
%     [pred_lab{p}, pred_colors(p,:), pred_styles(p)] = fn_get_label_styles(predictors{p},1);
% end
pred_colors = [178,24,43;...
              33,102,172;...
              244,165,130;...
              146,197,222]./256;
pred_styles  = {'-','-','--','--'};

%% Analysis Set Up
% Plotting parameters
plot_boxes = 1;
pow_y_lim = 1;%0.21;
font_size = 20;
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';

% Load SBJ info, stat info:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if ~isfield(st,'time_resolved') || ~st.time_resolved; error('only use time resovled stat_ids!'); end
if st.use_simon_tfr~=0; error('why use Simon TFR?'); end
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end
eval(['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m']);

% Load power band info
[theta_cf, betalo_cf, betahi_cf] = fn_get_sbj_peak_frequencies(SBJs,an_id);
theta_lim  = fn_compute_freq_lim(SBJs,theta_cf,'theta');
betalo_lim = fn_compute_freq_lim(SBJs,betalo_cf,'betalo');
betahi_lim = fn_compute_freq_lim(SBJs,betahi_cf,'betahi');

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

%% Load LMM results
if any(contains(predictors,'SV'))
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMM_timeresolved_SV.mat'];
    fig_suffix = '_SV';
else
    lmm_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_LMM_timeresolved.mat'];
    fig_suffix = '';
end
fprintf('Loading %s\n',lmm_fname);
load(lmm_fname);
if length(pow_vars)~=4; error('why not use PFC/BG and theta/beta?'); end
if (length(predictors)==4 && ~all(strcmp(predictors,{'reward_cur','effortS_cur','reward_prv','effortS_prv'}))) || ...
    (length(predictors)==2 && ~all(strcmp(predictors,{'SV_cur','SV_prv'})))
    error('predictor should be reward-effort or SV model!');
end

%% Extract plotting LME data
lmm_coef = nan([length(pow_vars) length(predictors) numel(plt_time_vec)]);
lmm_ci   = nan([length(pow_vars) length(predictors) 2 numel(plt_time_vec)]);
lmm_pval = nan([length(pow_vars) length(predictors) numel(plt_time_vec)]);
lmm_sig  = nan([length(pow_vars) length(predictors) numel(plt_time_vec)]);
for nrl_ix = 1:length(pow_vars)
    if ~all(strcmp(lmm_ts{nrl_ix,1}.CoefficientNames(2:end),predictors))
        error('check fixed effect coefficeint and order');
    end
    if ~strcmp(lmm_ts{nrl_ix,1}.ResponseName,pow_vars{nrl_ix})
        error('check response variable');
    end
    for p_ix = 1:length(predictors)
        coef_ix = strcmp(lmm_ts{nrl_ix,1}.CoefficientNames,predictors{p_ix});
        for t_ix = 1:length(plt_time_vec)
            lmm_coef(nrl_ix,p_ix,t_ix) = lmm_ts{nrl_ix,t_ix}.Coefficients.Estimate(coef_ix);
            lmm_pval(nrl_ix,p_ix,t_ix) = lmm_stat_ts{nrl_ix,t_ix,p_ix}.pValue(2);
            lmm_ci(nrl_ix,p_ix,1,t_ix) = lmm_ts{nrl_ix,t_ix}.Coefficients.Upper(coef_ix);%-lmm_coef(m_ix,t_ix);
            lmm_ci(nrl_ix,p_ix,2,t_ix) = lmm_ts{nrl_ix,t_ix}.Coefficients.Lower(coef_ix);%+lmm_coef(m_ix,t_ix);
        end
        
        % Correct for time points
        %     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
        [lmm_sig(nrl_ix,p_ix,:), ~, ~, ~]  = fdr_bh(lmm_pval(nrl_ix,p_ix,:));
    end
end

%% Get frequency and time ticks
time_ticks = [0:0.5:2];%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
time_tick_ix = nan(size(time_ticks));
time_tick_lab = cell(size(time_ticks));
for t = 1:numel(time_ticks)
    [~,time_tick_ix(t)] = min(abs(plt_time_vec-time_ticks(t)));
    time_tick_lab{t} = ['' num2str(time_ticks(t))];
end
if ~any(time_ticks==0); error('cant plot event with no tick at 0'); end

%% Plot coefficient time series
fig_name = ['GRP_LMM_ts_' an_id '_' stat_id fig_suffix];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for nrl_ix = 1:length(pow_vars)
    subplot(2,2,nrl_ix); hold on;
    line([plt_time_vec(1) plt_time_vec(end)],[0 0],'Color','k','LineWidth',1);
    % Plot main effect
    pred_lines = gobjects(size(predictors));
    for p_ix = 1:length(predictors)
        %     lmm_line = shadedErrorBar(plt_time_vec, squeeze(lmm_coef(m_ix,:)), ...
        %         squeeze(lmm_ci(m_ix,:,:)),'lineprops',{'Color',lmm_colors(m_ix,:),'LineWidth',1});
        pred_lines(p_ix) = plot(plt_time_vec,squeeze(lmm_coef(nrl_ix,p_ix,:)),...
            'Color',pred_colors(p_ix,:),'LineStyle',pred_styles{p_ix},'LineWidth',2);
%         p = patch([plt_time_vec flip(plt_time_vec)],[squeeze(lmm_ci(nrl_ix,p_ix,1,:))' flip(squeeze(lmm_ci(nrl_ix,p_ix,2,:)))'],...
%             lmm_colors(p_ix,:),'FaceAlpha',0.2);
%         p.EdgeColor = lmm_colors(p_ix,:);

        % Plot significant time periods
        if any(squeeze(lmm_sig(nrl_ix,p_ix,:)))
            % Find significant periods
            sig_chunks = fn_find_chunks(squeeze(lmm_sig(nrl_ix,p_ix,:)));
            sig_chunks(lmm_sig(nrl_ix,p_ix,sig_chunks(:,1))==0,:) = [];
            fprintf('%s %s -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                pow_vars{nrl_ix},predictors{p_ix},size(sig_chunks,1));

            % Plot Significance
            for sig_ix = 1:size(sig_chunks,1)
                if diff(sig_chunks(sig_ix,:))==0    % single windows
                    scatter(plt_time_vec(sig_chunks(sig_ix,1)),squeeze(lmm_coef(nrl_ix,p_ix,sig_chunks(sig_ix,1))),...
                        sig_scat_size,pred_colors(p_ix,:),'o','filled');
                else
                    line(plt_time_vec(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                        squeeze(lmm_coef(nrl_ix,p_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                        'Color',pred_colors(p_ix,:),'LineStyle',plt.sig_style,...
                        'LineWidth',plt.sig_width);
                end
            end
        end
    end
    
    ylims = ylim;
    if plot_boxes
        patch([0.5 0.5 1.5 1.5],[-pow_y_lim pow_y_lim pow_y_lim -pow_y_lim],'k','FaceAlpha',0,...
            'LineWidth',2,'LineStyle','--');
        box_str = '_boxes';
    else
        box_str = '';
    end
    
    % Axes and parameters
    title(pow_vars{nrl_ix}, 'interpreter', 'none');
    set(gca,'XLim', st.stat_lim);
    set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
    set(gca,'YLim',[-pow_y_lim pow_y_lim]);
    xlabel('Time (s)');
    ylabel('Coefficient');
    set(gca,'FontSize',font_size);
    if nrl_ix==1
        leg = legend(pred_lines,predictors,'location','best','interpreter','none','fontsize',14);
    end
end

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name box_str '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot coefficient time series for PFC theta ~ rew_prv and BG beta ~ effort_cur without error bars
fig_name = ['GRP_LMM_ts_' an_id '_' stat_id '_PFCtpR_BGbcE_noerrbr' fig_suffix];
figure('Name',fig_name,'units','normalized');%,...
%     'outerposition',[0 0 0.3 0.5],'Visible',fig_vis);
if contains(stat_id,'nrl0_out3main')
    pow_y_lim = [-0.4 0.7];
elseif contains(stat_id,'nrlz_out4main')
    pow_y_lim = [-0.15 0.15];
else, error('why not use the main 2 stat_ids?');
end    
[~,pR_colors,~] = fn_get_label_styles('reward_prv',3);
[~,cE_colors,~] = fn_get_label_styles('effortS_cur',3);

% Plot coefficient time series
hold on;
nrl_ix1 = strcmp(pow_vars,'PFC_theta'); nrl_ix2 = strcmp(pow_vars,'BG_betalo');
p_ix1 = strcmp(predictors,'reward_prv'); p_ix2 = strcmp(predictors,'effortS_cur');
line([0 2],[0 0],'Color','k','LineWidth',1);
t_line = plot(plt_time_vec,squeeze(lmm_coef(nrl_ix1,p_ix1,:)),'Color',pR_colors(2,:),'LineWidth',2);
b_line = plot(plt_time_vec,squeeze(lmm_coef(nrl_ix2,p_ix2,:)),'Color',cE_colors(2,:),'LineWidth',2);

ylims = ylim;
if plot_boxes
    patch([0.5 0.5 1.5 1.5],[pow_y_lim(1) pow_y_lim(2) pow_y_lim(2) pow_y_lim(1)],'k','FaceAlpha',0,...
        'LineWidth',2,'LineStyle','--');
    box_str = '_boxes';
else
    box_str = '';
end

% Axes and parameters
title('LMM effects', 'interpreter', 'none');
set(gca,'XLim', st.stat_lim);
set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
set(gca,'YLim',pow_y_lim);
xlabel('Time (s)');
ylabel('Coefficient');
legend([t_line,b_line],{'PFC Theta ~ Prev. Reward','BG Beta ~ Effort'},'Location','best');
set(gca,'FontSize',font_size);

% Save Figure
if save_fig
    fig_fname = [fig_dir fig_name box_str '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot coefficient time series by ROI without error bars
% fig_name = ['GRP_LMM_ts_' an_id '_' stat_id '_byROI_noerrbr' fig_suffix];
% figure('Name',fig_name,'units','normalized',...
%     'outerposition',[0 0 0.3 1],'Visible',fig_vis);
% pow_y_lim = 0.15;
% 
% % Plot PFC
% subplot(2,1,1); hold on;
% nrl_ix = strcmp(pow_vars,'PFC_theta');
% p_ix1 = strcmp(predictors,'reward_prv'); p_ix2 = strcmp(predictors,'effortS_cur');
% line([0 2],[0 0],'Color','k','LineWidth',1);
% t_line = plot(plt_time_vec,squeeze(lmm_coef(1,p_ix,:)),'Color',pred_colors(1,:),'LineWidth',2);
% b_line = plot(plt_time_vec,squeeze(lmm_coef(2,:)),'Color',pred_colors(2,:),'LineWidth',2);
% 
% ylims = ylim;
% if plot_boxes
%     patch([0.5 0.5 1.5 1.5],[-pow_y_lim pow_y_lim pow_y_lim -pow_y_lim],'k','FaceAlpha',0,...
%         'LineWidth',2,'LineStyle','--');
%     box_str = '_boxes';
% else
%     box_str = '';
% end
% 
% % Axes and parameters
% title('PFC LMM effects', 'interpreter', 'none');
% set(gca,'XLim', st.stat_lim);
% set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
% set(gca,'YLim',[-pow_y_lim pow_y_lim]);
% xlabel('Time (s)');
% ylabel('Coefficient');
% legend([t_line,b_line],{'Theta ~ Prev. Reward','Beta ~ Effort'},'Location','best');
% set(gca,'FontSize',font_size);
% 
% % Plot BG
% subplot(2,1,2); hold on;
% line([0 2],[0 0],'Color','k','LineWidth',1);
% t_line = plot(plt_time_vec,squeeze(lmm_coef(3,:)),'Color',pred_colors(3,:),'LineWidth',2);
% b_line = plot(plt_time_vec,squeeze(lmm_coef(4,:)),'Color',pred_colors(4,:),'LineWidth',2);
% 
% ylims = ylim;
% if plot_boxes
%     patch([0.5 0.5 1.5 1.5],[-pow_y_lim pow_y_lim pow_y_lim -pow_y_lim],'k','FaceAlpha',0,...
%         'LineWidth',2,'LineStyle','--');
%     box_str = '_boxes';
% else
%     box_str = '';
% end
% 
% % Axes and parameters
% title('BG LMM effects', 'interpreter', 'none');
% set(gca,'XLim', st.stat_lim);
% set(gca,'XTick', time_ticks);%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
% set(gca,'YLim',[-pow_y_lim pow_y_lim]);
% xlabel('Time (s)');
% ylabel('Coefficient');
% legend([t_line,b_line],{'Theta ~ Prev. Reward','Beta ~ Effort'},'Location','best');
% set(gca,'FontSize',font_size);
% 
% % Save Figure
% if save_fig
%     fig_fname = [fig_dir fig_name box_str '.' fig_ftype];
%     fprintf('Saving %s\n',fig_fname);
%     saveas(gcf,fig_fname);
% end
% 
