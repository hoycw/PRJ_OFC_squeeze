function fn_plot_LMM_gratton_bar(tbl,xvar_prv,xvar_cur,yvar)

% if ~contains(xvar_prv,'prv') || ~contains(xvar_cur,'cur')
%     error('input xvar must have cur/prv for gratton violin plot');
% end

if contains(xvar_prv,'_prv')
    prv_lab = ['Previous ' strrep(xvar_prv,'_prv','')];
elseif contains(xvar_prv,'_cur')
    prv_lab = ['Current ' strrep(xvar_prv,'_cur','')];
end
if contains(xvar_cur,'_prv')
    cur_lab = ['Previous ' strrep(xvar_cur,'_prv','')];
elseif contains(xvar_cur,'_cur')
    cur_lab = ['Current ' strrep(xvar_cur,'_cur','')];
end
lo_color = [1 1 1];%[253,184,99]./256;
hi_color = [0.7 0.7 0.7];%[178,171,210]./256;
scat_sz  = 50;

if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end

%% Compute summary stats for combinations of previous and current
if any(contains(xvar_prv,{'reward','effort'}))
    prv_cut_off = [-0.5 0.5];
else
    prv_cut_off = [median(tbl.(xvar_prv)) median(tbl.(xvar_prv))];
end
prv_lo_idx = tbl.(xvar_prv)<prv_cut_off(1);
prv_hi_idx = tbl.(xvar_prv)>prv_cut_off(2);
if any(contains(xvar_cur,{'reward','effort'}))
    cur_cut_off = [-0.5 0.5];
else
    cur_cut_off = [median(tbl.(xvar_cur)) median(tbl.(xvar_cur))];
end
cur_lo_idx = tbl.(xvar_cur)<cur_cut_off(1);
cur_hi_idx = tbl.(xvar_cur)>cur_cut_off(2);
cond_labs = {'lPlC','lPhC','hPlC','hPhC'};
cond_idx(:,1) = prv_lo_idx & cur_lo_idx;
cond_idx(:,2) = prv_lo_idx & cur_hi_idx;
cond_idx(:,3) = prv_hi_idx & cur_lo_idx;
cond_idx(:,4) = prv_hi_idx & cur_hi_idx;

grp_means = nan(size(cond_labs));
grp_sems  = nan(size(cond_labs));
for c = 1:length(cond_labs)
    grp_means(c) = mean(tbl.(yvar)(cond_idx(:,c)));
    grp_sems(c)  = std(tbl.(yvar)(cond_idx(:,c)))./sqrt(sum(cond_idx(:,c)));
end
% Print means:
fprintf('Difference in %s for %s:%s:\n',yvar,xvar_cur,xvar_prv);
fprintf('\tLow %s: %.2f high %s - %.2f low %s = %.3f\n',xvar_prv,grp_means(2),...
    xvar_cur,grp_means(1),xvar_cur,grp_means(2)-grp_means(1));
fprintf('\tHi  %s: %.2f high %s - %.2f low %s = %.3f\n',xvar_prv,grp_means(4),...
    xvar_cur,grp_means(3),xvar_cur,grp_means(4)-grp_means(3));

%% Plot error bars
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_prv '-' xvar_cur '_gratton_bar'];
figure('Name',fig_name); hold on;

% Plot group means as bars
cond_xpos = [0.5 1.5 3.5 4.5];
bars = bar(cond_xpos,diag(grp_means),'stacked');
for c = 1:4
    if contains(cond_labs{c},'lC'); plt_color = lo_color; else; plt_color = hi_color; end
    set(bars(c),'FaceColor',plt_color,'EdgeColor','k');
    errorbar(cond_xpos(c),grp_means(c),grp_sems(c),...
        'Color','k','LineWidth',3);
end
ylabel(yvar);
if contains(yvar,'decision')
    ylim([0 1]); yticks(0:0.2:1);
else
    ylim([min(grp_means(:))-range(grp_means(:))*0.1 max(grp_means(:))+range(grp_means(:))*0.1]);
end
set(gca,'XTick',[1 4]);
set(gca,'XTickLabel',{['Low ' prv_lab],['High ' prv_lab]});
title(['Effects of ' prv_lab ' x ' cur_lab]);
set(gca,'FontSize',16);
legend(bars(1:2),{['Low ' cur_lab],['High ' cur_lab]},'Location','best');

end