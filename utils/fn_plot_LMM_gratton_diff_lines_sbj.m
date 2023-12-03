function fn_plot_LMM_gratton_diff_lines_sbj(tbl,xvar_prv,xvar_cur,yvar)

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
plt_color = [0.7 0.7 0.7];
scat_sz  = 50;

if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

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

% Compute means within condition
grp_means = nan(size(cond_labs));
grp_sems  = nan(size(cond_labs));
sbj_means = nan([4 length(cond_labs)]);
for c = 1:length(cond_labs)
    grp_means(c) = mean(tbl.(yvar)(cond_idx(:,c)));
    grp_sems(c)  = std(tbl.(yvar)(cond_idx(:,c)))./sqrt(sum(cond_idx(:,c)));
    for s = 1:4
        sbj_means(s,c) = mean(tbl.(yvar)(cond_idx(:,c) & tbl.sbj_n==s));
    end
end

% Compute differences between conditions
sbj_diff  = nan([4 2]);
for s = 1:4
    sbj_diff(s,1) = sbj_means(s,2)-sbj_means(s,1);
    sbj_diff(s,2) = sbj_means(s,4)-sbj_means(s,3);
end
grp_diff(1) = mean(sbj_diff(:,1));
grp_diff_sem(1) = std(sbj_diff(:,1))./sqrt(4);
grp_diff(2) = mean(sbj_diff(:,2));
grp_diff_sem(2) = std(sbj_diff(:,2))./sqrt(4);

% Print means:
fprintf('Difference in %s for %s:%s:\n',yvar,xvar_cur,xvar_prv);
fprintf('\tLow %s: %.2f high %s - %.2f low %s = %.3f\n',xvar_prv,grp_means(2),...
    xvar_cur,grp_means(1),xvar_cur,grp_diff(1));
fprintf('\tHi  %s: %.2f high %s - %.2f low %s = %.3f\n',xvar_prv,grp_means(4),...
    xvar_cur,grp_means(3),xvar_cur,grp_diff(2));

%% Plot error bars
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_prv '-' xvar_cur '_gratton_diff_lines_sbj'];
figure('Name',fig_name,'units','norm','outerposition',[0 0 0.2 0.5]); hold on;

% Plot group means as bars
cond_xpos = [0.5 1.5];
sbj_xvals = linspace(-0.05,0.05,4);
bars = bar(cond_xpos,diag(grp_diff),'stacked');
for c = 1:2
    set(bars(c),'FaceColor',plt_color,'EdgeColor','k');
    errorbar(cond_xpos(c),grp_diff(c),grp_diff_sem(c),'Color','k','LineWidth',3);
end
% Plot subject mean differences as lines
for s = 1:4
    scatter(sbj_xvals(s)+cond_xpos(1),sbj_diff(s,1),scat_sz,sbj_colors(s,:));
    scatter(sbj_xvals(s)+cond_xpos(2),sbj_diff(s,2),scat_sz,sbj_colors(s,:));
    line(sbj_xvals(s)+cond_xpos,sbj_diff(s,:),'Color',sbj_colors(s,:),'LineWidth',2);
end
ylabel([yvar ' hi-lo ' xvar_prv]);
ylim([min(sbj_diff(:))-range(sbj_diff(:))*0.1 max(sbj_diff(:))+range(sbj_diff(:))*0.1]);
set(gca,'XTick',cond_xpos);
set(gca,'XTickLabel',{['Low ' prv_lab],['High ' prv_lab]});
title([cur_lab ' diffs for high/low ' prv_lab]);
set(gca,'FontSize',16);

end