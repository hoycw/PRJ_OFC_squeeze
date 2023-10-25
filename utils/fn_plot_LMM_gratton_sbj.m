function fn_plot_LMM_gratton_sbj(SBJs,tbl,xvar,yvar)
%% Gratton line plot with single subejct overlays

% Set up and check vaeriables
if contains(xvar,'prv') || contains(xvar,'cur')
    error('input xvar name without cur or prv for gratton plot');
end
if ~all(strcmp(SBJs,{'PFC03','PFC04','PFC05','PFC01'})); error('SBJs wrong'); end
if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end

sbj_xvals = [-0.1 -0.05 0.05 0.1];
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta

%% Compute summary stats and plot per SBJ
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar '_gratton_sbj'];
figure('Name',fig_name); hold on;

for s = 1:length(SBJs)
    % Compute summary stats for combinations of previous and current
    prv_low_idx = tbl.([xvar '_prv'])<median(tbl.([xvar '_prv'])) & tbl.sbj_n==s;
    cur_low_idx = tbl.([xvar '_cur'])<median(tbl.([xvar '_cur'])) & tbl.sbj_n==s;
    lL_avg = mean(tbl.(yvar)(prv_low_idx & cur_low_idx));
    lH_avg = mean(tbl.(yvar)(prv_low_idx & ~cur_low_idx));
    hL_avg = mean(tbl.(yvar)(~prv_low_idx & cur_low_idx));
    hH_avg = mean(tbl.(yvar)(~prv_low_idx & ~cur_low_idx));
    lL_sem = std(tbl.(yvar)(prv_low_idx & cur_low_idx))./sqrt(sum(prv_low_idx & cur_low_idx));
    lH_sem = std(tbl.(yvar)(prv_low_idx & ~cur_low_idx))./sqrt(sum(prv_low_idx & ~cur_low_idx));
    hL_sem = std(tbl.(yvar)(~prv_low_idx & cur_low_idx))./sqrt(sum(~prv_low_idx & cur_low_idx));
    hH_sem = std(tbl.(yvar)(~prv_low_idx & ~cur_low_idx))./sqrt(sum(~prv_low_idx & ~cur_low_idx));
    
    % Plot error bars
    cur_lo_line = errorbar([1 2]+sbj_xvals(s),[lL_avg hL_avg],[lL_sem hL_sem],'Color',sbj_colors(s,:),...
        'LineStyle','--','LineWidth',1.5);
    cur_hi_line = errorbar([1 2]+sbj_xvals(s),[lH_avg hH_avg],[lH_sem hH_sem],'Color',sbj_colors(s,:),...
        'LineStyle','-','LineWidth',1.5);
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

% Plot error bars
cur_lo_line = errorbar([1 2],[lL_avg hL_avg],[lL_sem hL_sem],'Color','k',...
    'LineStyle','--','LineWidth',3);
cur_hi_line = errorbar([1 2],[lH_avg hH_avg],[lH_sem hH_sem],'Color','k',...
    'LineStyle','-','LineWidth',3);
ylabel(yvar);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{['Low Previous ' xvar],['High Previous ' xvar]});
xlim([0.5 2.5]);
legend([cur_lo_line, cur_hi_line],{['Low Current ' xvar],['High Current ' xvar]},'Location','best');
title(['Effects of Previous/Current ' xvar]);
set(gca,'FontSize',16);

end