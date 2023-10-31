function fn_plot_LMM_interaction_partial_dependence(lme,xvar_plt,xvar_div)
% Line plots of yvar as a function of xvar_plt, discretized by xvar_div
% (one line per value of xvar_div)

%% Set up and check variables
[xplt_lab, ~, ~] = fn_get_label_styles(xvar_plt,1);
[xdiv_lab, ~, ~] = fn_get_label_styles(xvar_div,1);
[yvar_lab, ~, ~] = fn_get_label_styles(lme.ResponseName,1);

%% Compute Partial dependence
% div_qs = quantile(tbl.(xvar_div),n_quant_div);
% plt_qs = linspace(min(tbl.(xvar_plt)),max(tbl.(xvar_plt)),n_quant_plt);
% plt_pts = nan([n_quant_div*n_quant_plt 2]);
% i = 0;
% for xd = 1:n_quant_div
%     for xp = 1:n_quant_plt
%         i = i+1;
%         plt_pts(i,:) = [div_qs(xd) plt_qs(xp)];
%     end
% end
figure;
pd_plt = plotPartialDependence(lme,{xvar_div xvar_plt},'UseParallel',true);%,'QueryPoints',plt_pts);
view(0,90);
pd_zvals = pd_plt.Children.CData;
pd_xvals = pd_plt.Children.XData(1,:);
pd_yvals = pd_plt.Children.YData(:,1);
close(gcf);

%% Plot partial dependence curves
fig_name = ['GRP_TFR_LMM_results_' lme.ResponseName '_' xvar_plt ...
    '_by_' xvar_div '_partial_dependence_interaction_mat'];
figure('Name',fig_name); hold on;
imagesc(pd_xvals,pd_yvals,pd_zvals);

xlabel(xdiv_lab);
xlim([min(pd_xvals) max(pd_xvals)]);
ylabel(xplt_lab);
ylim([min(pd_yvals) max(pd_yvals)]);
cbar = colorbar;
cbar.Label.String = yvar_lab;
title(['Interaction of ' xplt_lab{1} ' vs ' strrep(xdiv_lab{1},' Q1','')]);
set(gca,'FontSize',16);

end