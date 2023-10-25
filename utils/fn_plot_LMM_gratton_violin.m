function fn_plot_LMM_gratton_violin(tbl,xvar_prv,xvar_cur,yvar)

if ~contains(xvar_prv,'prv') || ~contains(xvar_cur,'cur')
    error('input xvar must have cur/prv for gratton violin plot');
end
if contains(yvar,'decision'); error('violins of bianry variables make no sense'); end

prv_lab = strrep(xvar_prv,'_prv','');
cur_lab = strrep(xvar_cur,'_cur','');
lo_color = [230,97,1]./256;
hi_color = [94,60,153]./256;

%% Compute summary stats for combinations of previous and current
prv_low_idx = tbl.(xvar_prv)<median(tbl.(xvar_prv));
cur_low_idx = tbl.(xvar_cur)<median(tbl.(xvar_cur));
cond_labs = {'lPlC','lPhC','hPlC','hPhC'};
cond_idx(:,1) = prv_low_idx & cur_low_idx;
cond_idx(:,2) = prv_low_idx & ~cur_low_idx;
cond_idx(:,3) = ~prv_low_idx & cur_low_idx;
cond_idx(:,4) = ~prv_low_idx & ~cur_low_idx;

v_data = struct;
% sbj_means = nan([length(SBJs) length(cond_labs)]);
for c = 1:length(cond_labs)
    v_data.(cond_labs{c}) = tbl.(yvar)(cond_idx(:,c));
%     for s = 1:length(SBJs)
%         sbj_means(s,c) = mean(tbl.(yvar)(cond_idx(:,c) & strcmp(tbl.sbj_id,SBJs{s})));
%     end
end

%% Plot violins
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_prv '-' xvar_cur '_gratton_violin'];
figure('Name',fig_name); hold on;

violins = violinplot(v_data);
for p = 1:length(cond_labs)
    violins(p).MeanPlot.LineWidth = 3;
    violins(p).MeanPlot.Visible   = 'on';
    violins(p).ScatterPlot.Visible = 'off';
    if contains(cond_labs{p},'lC')
        violins(p).ViolinColor = lo_color;
    else
        violins(p).ViolinColor = hi_color;
    end
    % Violin mean line is off since ksdensity is tossing some of the
    % values for some reason, so manually set it and adjust width
    violins(p).MeanPlot.YData = repmat(mean(v_data.(cond_labs{p})),1,2);
    [~,mean_ix] = min(abs(violins(p).ViolinPlot.YData-mean(v_data.(cond_labs{p}))));
    violins(p).MeanPlot.XData = [p-(violins(p).ViolinPlot.XData(mean_ix)-p), ...
        violins(p).ViolinPlot.XData(mean_ix)];
    violins(p).MeanPlot.LineWidth = 4;
end
% % Add Subject means
% sbj_xvals = linspace(-0.1,0.1,length(SBJs));
% sbj_lines = gobjects(size(SBJs));
% for i = 1:2
%     for s = 1:length(SBJs)
%         sbj_lines(s) = line(([sbj_xvals(s)+1 sbj_xvals(s)+2]+(i-1)*2),...
%             [sbj_means(s,1+(i-1)*2) sbj_means(s,2+(i-1)*2)],'Color',sbj_colors(s,:),'LineWidth',2);
%     end
% end
ylabel(yvar);
xticks([1.5 3.5]);
xticklabels({['Low Previous ' prv_lab],['High Previous ' prv_lab]});
legend([violins(1).ViolinPlot;violins(2).ViolinPlot],...
    {['Low Current ' cur_lab],['High Current ' cur_lab]},'location','best');
title(['Effects of Previous ' prv_lab ':Current ' cur_lab]);
ax = gca;
ax.FontSize = 16;
ax.XLabel.Interpreter = 'none';

end