function fn_plot_LMM_gratton(tbl,xvar_prv,xvar_cur,yvar)

if ~contains(xvar_prv,'prv') || ~contains(xvar_cur,'cur')
    error('input xvar must have cur/prv for gratton violin plot');
end

prv_lab = strrep(xvar_prv,'_prv','');
cur_lab = strrep(xvar_cur,'_cur','');
lo_color = [230,97,1]./256;
hi_color = [94,60,153]./256;
scat_sz  = 40;

if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end

sbj_xvals = linspace(-0.05,0.05,4);
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

%% Compute summary stats for combinations of previous and current
prv_low_idx = tbl.(xvar_prv)<median(tbl.(xvar_prv));
cur_low_idx = tbl.(xvar_cur)<median(tbl.(xvar_cur));
cond_labs = {'lPlC','lPhC','hPlC','hPhC'};
cond_idx(:,1) = prv_low_idx & cur_low_idx;
cond_idx(:,2) = prv_low_idx & ~cur_low_idx;
cond_idx(:,3) = ~prv_low_idx & cur_low_idx;
cond_idx(:,4) = ~prv_low_idx & ~cur_low_idx;

% lL_avg = mean(tbl.(yvar)(prv_low_idx & cur_low_idx));
% lH_avg = mean(tbl.(yvar)(prv_low_idx & ~cur_low_idx));
% hL_avg = mean(tbl.(yvar)(~prv_low_idx & cur_low_idx));
% hH_avg = mean(tbl.(yvar)(~prv_low_idx & ~cur_low_idx));
% lL_sem = std(tbl.(yvar)(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
% lH_sem = std(tbl.(yvar)(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
% hL_sem = std(tbl.(yvar)(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
% hH_sem = std(tbl.(yvar)(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));

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

%% Plot error bars
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_prv '-' xvar_cur '_gratton'];
figure('Name',fig_name); hold on;

for c = 1:4
    if contains(cond_labs{c},'lC'); scat_color = lo_color; else; scat_color = hi_color; end
    if contains(cond_labs{c},'lP'); xpos = 1; else; xpos = 2; end
    scatter(sbj_xvals+xpos,sbj_means(:,c),scat_sz,scat_color);
end
cur_lo_line = errorbar([1 2],grp_means(contains(cond_labs,'lC')),grp_sems(contains(cond_labs,'lP')),...
    'Color',lo_color,'LineWidth',2);
cur_hi_line = errorbar([1 2],grp_means(contains(cond_labs,'hC')),grp_sems(contains(cond_labs,'hP')),...
    'Color',hi_color,'LineWidth',2);
ylabel(yvar);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{['Low Previous ' prv_lab],['High Previous ' prv_lab]});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{['Low Current ' cur_lab],['High Current ' cur_lab]},'Location','best');
title(['Effects of Previous ' prv_lab ':Current ' cur_lab]);
set(gca,'FontSize',16);

end