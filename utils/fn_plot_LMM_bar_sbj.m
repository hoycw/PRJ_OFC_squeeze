function fn_plot_LMM_bar_sbj(tbl,xvar,yvar)

% if ~contains(xvar_prv,'prv') || ~contains(xvar_cur,'cur')
%     error('input xvar must have cur/prv for gratton violin plot');
% end

if contains(xvar,'_prv')
    x_lab = ['Previous ' strrep(xvar,'_prv','')];
elseif contains(xvar,'_cur')
    x_lab = ['Current ' strrep(xvar,'_cur','')];
end
lo_color = [1 1 1];%[253,184,99]./256;
hi_color = [0.7 0.7 0.7];%[178,171,210]./256;
scat_sz  = 50;

if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

%% Compute summary stats for levels of xvar
if any(contains(xvar,{'reward','effort'}))
    cut_off = [-0.5 0.5];
else
    cut_off = [median(tbl.(xvar)) median(tbl.(xvar))];
end
lo_idx = tbl.(xvar)<cut_off(1);
hi_idx = tbl.(xvar)>cut_off(2);
cond_labs = {'lo','hi'};
cond_idx(:,1) = lo_idx;
cond_idx(:,2) = hi_idx;

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
% Print means:
fprintf('Difference in %s for %s:\n',yvar,xvar);
fprintf('\t%s: high %.2f - low %.2f = %.3f\n',xvar,grp_means(2),...
    grp_means(1),grp_means(2)-grp_means(1));

%% Plot error bars
% fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar '_bar_sbj'];
% figure('Name',fig_name); hold on;

% Plot group means as bars
cond_xpos = [1 2];
sbj_xvals = linspace(-0.05,0.05,4);
bars = bar(cond_xpos,diag(grp_means),'stacked');
for s = 1:4
    for c= 1:2
        scatter(sbj_xvals(s)+cond_xpos(c),sbj_means(s,c),scat_sz,sbj_colors(s,:));
    end
    line(sbj_xvals(s)+cond_xpos(1:2),sbj_means(s,1:2),'Color',sbj_colors(s,:),'LineWidth',2);
end
for c = 1:2
    if contains(cond_labs{c},'lo'); plt_color = lo_color; else; plt_color = hi_color; end
    set(bars(c),'FaceColor',plt_color,'EdgeColor','k');
    errorbar(cond_xpos(c),grp_means(c),grp_sems(c),...
        'Color','k','LineWidth',3);
end
ylabel(yvar);
if contains(yvar,'decision')
    ylim([0 1]); yticks(0:0.2:1);
else
    ylim([min(sbj_means(:))-range(sbj_means(:))*0.1 max(sbj_means(:))+range(sbj_means(:))*0.1]);
end
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{['Low ' x_lab],['High ' x_lab]});
xlim([0 3]);
title(['Effects of ' x_lab]);
set(gca,'FontSize',16);

end