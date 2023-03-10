%% Plot group model coeffifient TFR matrix per circular-linear regressor + R2
%   Betas are outlined by significance
%   R2 is plotted as separatefigure, one matrix per regressor

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');

close all
clear all

%% Analysis parameters:
an_id = 'TFRmth_S03t2_f2t30_fourier'; stat_id = 'S5t15_bhvz_nrlphs';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

% Plotting parameters
% symmetric_clim = 1;
font_size = 18;
save_fig  = 1;
fig_ftype = 'png';
fig_vis   = 'on';
x_step    = 0.2;
ns_alpha  = 0.4;

%% Prep stuff
if contains(an_id,'f2t30'); freq_ticks = 5:5:30; else; freq_ticks = 5:5:35; end
% if symmetric_clim; clim_str = '_sym'; else clim_str = ''; end
if contains(an_id,'_S') || contains(an_id,'simon')
    plt_id = 'ts_S2t2_evnts_sigLine';
elseif contains(an_id,'_D')
    plt_id = 'ts_D1t1_evnts_sigLine';
else
    error('couldnt pick plt_id based on an_id');
end

fig_dir   = [prj_dir 'results/TFR/' an_id '/CLreg_' stat_id '/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

plt_vars_cmd = ['run ' prj_dir 'scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
if contains(an_id,'_S')
    time_ticks = [0:0.5:2];%plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
elseif contains(an_id,'_D')
    time_ticks = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
end

%% Load data
results_fname = [prj_dir 'data/GRP/GRP_CLreg_' an_id '_' stat_id '.mat'];
fprintf('\tLoading %s...\n',results_fname);
load(results_fname,'betas','r2s','qvals','pvals','reg_vars','time_vec','freq_vec');

%% Plot results with FDR correction
% Get significance mask
sig_mask = ones(size(qvals))*ns_alpha;
sig_mask(qvals<=0.05) = 1;

for reg_ix = 1:size(reg_vars,1)
    % Get color limits across all regressors
    vals = betas(reg_ix,:,:);
    r2vals = r2s(reg_ix,:,:);
    beta_clim = [-max(abs(vals(:))) max(abs(vals(:)))];
    r2_clim   = [min(r2vals(:)) max(r2vals(:))];
    
    %% Create Beta Plot
    fig_name = ['GRP_CLreg_' an_id '_' stat_id '_' reg_vars{reg_ix,1} '_' reg_vars{reg_ix,2} '_qmask'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Plot coefficient matrix
    subplot(1,2,1); hold on;
    im = imagesc(time_vec, freq_vec, squeeze(betas(reg_ix,:,:)), beta_clim);
    im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
    
    line([0 0],ylim,'Color','k');
    
    % Axes and parameters
    set(gca,'YDir','normal');
    set(gca,'YLim',[min(freq_vec) max(freq_vec)]);
    set(gca,'XLim',[min(time_vec) max(time_vec)]);
    set(gca,'XTick',[min(time_vec):x_step:max(time_vec)]);
    title([reg_vars{reg_ix,1} ': Phase-' reg_vars{reg_ix,2} ' Beta'],'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',font_size);
    
    % Model Fit R2 Plots
    subplot(1,2,2); hold on;
    im = imagesc(time_vec, freq_vec, squeeze(r2s(reg_ix,:,:)),r2_clim);
    im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
    
    line([0 0],ylim,'Color','k');
    
    set(gca,'YDir','normal');
    set(gca,'YLim',[min(freq_vec) max(freq_vec)]);
    set(gca,'XLim',[min(time_vec) max(time_vec)]);
    set(gca,'XTick',[min(time_vec):x_step:max(time_vec)]);
    title([reg_vars{reg_ix,1} ': Phase-' reg_vars{reg_ix,2} ' R2'],'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',font_size);
    
%     %% Report max beta stats
%     for reg_ix = 1:numel(reg_lab)
%         max_beta = 0; max_f_ix = 0; max_t_ix = 0;
%         for f_ix = 1:numel(freq_vec)
%             if max(abs(betas(reg_ix,f_ix,:))) > abs(max_beta)
%                 max_tmp = max(betas(reg_ix,f_ix,:));
%                 min_tmp = min(betas(reg_ix,f_ix,:));
%                 if abs(max_tmp) > abs(min_tmp)
%                     [max_beta, max_t_ix] = max(betas(reg_ix,f_ix,:));
%                 else
%                     [max_beta, max_t_ix] = min(betas(reg_ix,f_ix,:));
%                 end
%                 max_f_ix = f_ix;
%             end
%         end
%         fprintf('%s max beta = %.03f at %.03f s and %.03f Hz; p = %.20f\n',reg_lab{reg_ix},max_beta,...
%             time_vec(max_t_ix),freq_vec(max_f_ix),qvals(reg_ix,max_f_ix,max_t_ix));
%     end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

%% Plot results with no correction for multiple comparisons
% Get significance mask
sig_mask = ones(size(pvals))*ns_alpha;
sig_mask(pvals<=0.05) = 1;

for reg_ix = 1:size(reg_vars,1)
    % Get color limits across all regressors
    vals = betas(reg_ix,:,:);
    r2vals = r2s(reg_ix,:,:);
    beta_clim = [-max(abs(vals(:))) max(abs(vals(:)))];
    r2_clim   = [min(r2vals(:)) max(r2vals(:))];
    
    %% Create Beta Plot
    fig_name = ['GRP_CLreg_' an_id '_' stat_id '_' reg_vars{reg_ix,1} '_' reg_vars{reg_ix,2} '_pmask'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Plot coefficient matrix
    subplot(1,2,1); hold on;
    im = imagesc(time_vec, freq_vec, squeeze(betas(reg_ix,:,:)), beta_clim);
    im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
    
    line([0 0],ylim,'Color','k');
    
    % Axes and parameters
    set(gca,'YDir','normal');
    set(gca,'YLim',[min(freq_vec) max(freq_vec)]);
    set(gca,'XLim',[min(time_vec) max(time_vec)]);
    set(gca,'XTick',[min(time_vec):x_step:max(time_vec)]);
    title([reg_vars{reg_ix,1} ': Phase-' reg_vars{reg_ix,2} ' Beta'],'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',font_size);
    
    % Model Fit R2 Plots
    subplot(1,2,2); hold on;
    im = imagesc(time_vec, freq_vec, squeeze(r2s(reg_ix,:,:)),r2_clim);
    im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
    
    line([0 0],ylim,'Color','k');
    
    set(gca,'YDir','normal');
    set(gca,'YLim',[min(freq_vec) max(freq_vec)]);
    set(gca,'XLim',[min(time_vec) max(time_vec)]);
    set(gca,'XTick',[min(time_vec):x_step:max(time_vec)]);
    title([reg_vars{reg_ix,1} ': Phase-' reg_vars{reg_ix,2} ' R2'],'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',font_size);
    
%     %% Report max beta stats
%     for reg_ix = 1:numel(reg_lab)
%         max_beta = 0; max_f_ix = 0; max_t_ix = 0;
%         for f_ix = 1:numel(freq_vec)
%             if max(abs(betas(reg_ix,f_ix,:))) > abs(max_beta)
%                 max_tmp = max(betas(reg_ix,f_ix,:));
%                 min_tmp = min(betas(reg_ix,f_ix,:));
%                 if abs(max_tmp) > abs(min_tmp)
%                     [max_beta, max_t_ix] = max(betas(reg_ix,f_ix,:));
%                 else
%                     [max_beta, max_t_ix] = min(betas(reg_ix,f_ix,:));
%                 end
%                 max_f_ix = f_ix;
%             end
%         end
%         fprintf('%s max beta = %.03f at %.03f s and %.03f Hz; p = %.20f\n',reg_lab{reg_ix},max_beta,...
%             time_vec(max_t_ix),freq_vec(max_f_ix),qvals(reg_ix,max_f_ix,max_t_ix));
%     end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end
