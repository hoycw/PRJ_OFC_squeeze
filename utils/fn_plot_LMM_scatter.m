function fn_plot_LMM_scatter(SBJs,tbl,xvar,yvar,lmm,pval)

x_fudge = 0.2;
scat_sz = 20;
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta
coef_ix = strcmp(lmm.CoefficientNames,xvar);

fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar '_scatter'];
figure('Name',fig_name); hold on;
xvals = min(tbl.(xvar))-x_fudge:0.01:max(tbl.(xvar))+x_fudge;
for s = 1:length(SBJs)
    scatter(tbl.(xvar)(tbl.sbj_n==s),tbl.(yvar)(tbl.sbj_n==s),...
        scat_sz,'k');%sbj_colors(s,:));
    hold on;
    mdl = fitlm(tbl.(xvar)(tbl.sbj_n==s),tbl.(yvar)(tbl.sbj_n==s));
    yvals = mdl.Coefficients.Estimate(1) + xvals*mdl.Coefficients.Estimate(2);
    line(xvals,yvals,'Color',sbj_colors(s,:),'LineWidth',3,'LineStyle',':');
end
yvals = lmm.Coefficients.Estimate(1) + xvals*lmm.Coefficients.Estimate(coef_ix);
line(xvals,yvals,'Color','k','LineWidth',5);
xlabel(xvar);
ylabel(yvar);
title(['LMM coef = ' num2str(lmm.Coefficients.Estimate(coef_ix)) '; p = ' num2str(pval,'%.03f')]);
set(gca,'FontSize',16);

end