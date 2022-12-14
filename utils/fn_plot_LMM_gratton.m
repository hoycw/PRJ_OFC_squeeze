function fn_plot_LMM_gratton(tbl,xvar,yvar)

if contains(xvar,'prv') || contains(xvar,'cur')
    error('input xvar name without cur or prv for gratton plot');
end

%% Compute summary stats for combinations of previous and current
prv_low_idx = tbl.([xvar '_prv'])<median(tbl.([xvar '_prv']));
cur_low_idx = tbl.([xvar '_cur'])<median(tbl.([xvar '_cur']));
lL_avg = mean(tbl.(yvar)(prv_low_idx & cur_low_idx));
lH_avg = mean(tbl.(yvar)(prv_low_idx & ~cur_low_idx));
hL_avg = mean(tbl.(yvar)(~prv_low_idx & cur_low_idx));
hH_avg = mean(tbl.(yvar)(~prv_low_idx & ~cur_low_idx));
lL_sem = std(tbl.(yvar)(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
lH_sem = std(tbl.(yvar)(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
hL_sem = std(tbl.(yvar)(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
hH_sem = std(tbl.(yvar)(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));

%% Plot error bars
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar '_gratton'];
figure('Name',fig_name); hold on;
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','b','LineWidth',2);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','r','LineWidth',2);
ylabel(yvar);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{['Low Previous ' xvar],['High Previous ' xvar]});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{['Low Current ' xvar],['High Current ' xvar]},'Location','best');
title(['Effects of Previous/Current ' xvar]);
set(gca,'FontSize',16);

end