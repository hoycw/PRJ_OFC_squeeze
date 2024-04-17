function fn_plot_LMM_quantile_line_interaction(tbl,xvar_plt,xvar_div,yvar,n_quant_div,n_quant_plt,varargin)
% Line plots of yvar as a function of xvar_plt, discretized by xvar_div
% (one line per value of xvar_div)

%% Set up and check variables
scat_sz  = 50;
x_fudge = 0.2;

sbj_ns = unique(tbl.sbj_n);
% if ~all(sbj_ns'==[1 2 3 4]); error('SBJs in tbl mismatch'); end
[xplt_lab, ~, ~] = fn_get_label_styles(xvar_plt,1);
[xdiv_lab, xdiv_colors, xdiv_styles] = fn_get_label_styles(xvar_div,n_quant_div);
[yvar_lab, ~, ~] = fn_get_label_styles(yvar,1);

%% Discretize predictors into quantiles
xdiv_idx = zeros(size(tbl.sbj_n));
xplt_idx = zeros(size(tbl.sbj_n));
for s_ix = 1:length(sbj_ns)
    s = sbj_ns(s_ix);
    % Get bin edges
%     if continuous_data
        div_qs = quantile(tbl.(xvar_div)(tbl.sbj_n==s),n_quant_div);
        plt_qs = quantile(tbl.(xvar_plt)(tbl.sbj_n==s),n_quant_plt);
%     else
%         div_qs = unique(tbl.(xvar_div)(tbl.sbj_n==s));
%         plt_qs = unique(tbl.(xvar_plt)(tbl.sbj_n==s));
%         if length(div_qs)~=n_quant_div
%             warning(['requested ' num2str(n_quant_div) ' quantiles but ' num2str(length(div_qs)) ' are in the data']);
%         end
%         if length(plt_qs)~=n_quant_plt
%             warning(['requested ' num2str(n_quant_plt) ' quantiles but ' num2str(length(plt_qs)) ' are in the data']);
%         end
%     end
    div_edges = nan([n_quant_div-1 1]);
    for i = 1:n_quant_div-1
        div_edges(i) = mean(div_qs(i:i+1));
    end
    plt_edges = nan([n_quant_plt-1 1]);
    for i = 1:n_quant_plt-1
        plt_edges(i) = mean(plt_qs(i:i+1));
    end
    
    % Assign to quantiles
    if n_quant_div==2 && any(contains(xvar_div,{'reward','effort'}))
        % Override median split for reward/effort to be high/low split
        if length(unique(tbl.(xvar_div)))==5 && all(unique(tbl.(xvar_div))'==[1 4 7 10 13])
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<5) = 1;
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)>8)  = 2;
        elseif length(unique(tbl.(xvar_div)))==5 && all(unique(tbl.(xvar_div))'==[0.16 0.32 0.48 0.64 0.8])
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<0.4) = 1;
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)>0.5)  = 2;
        else
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<-0.5) = 1;
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)>0.5)  = 2;
        end
    else
        for i = 1:n_quant_div
            if i==1                 % first quantile
                xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<div_edges(i)) = i;
            elseif i==n_quant_div   % last quantile
                xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)>=div_edges(i-1)) = i;
            else                    % middle quantiles
                xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<div_edges(i) &...
                    tbl.(xvar_div)>=div_edges(i-1)) = i;
            end
        end
    end
    for i = 1:n_quant_plt
        if i==1                 % first quantile
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)<plt_edges(i)) = i;
        elseif i==n_quant_plt   % last quantile
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)>=plt_edges(i-1)) = i;
        else                    % middle quantiles
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)<plt_edges(i) &...
                tbl.(xvar_plt)>=plt_edges(i-1)) = i;
        end
    end
end

%% Average within bins
x_mn = nan([n_quant_plt n_quant_div]);
y_mn = nan([n_quant_plt n_quant_div]);
y_se = nan([n_quant_plt n_quant_div]);
for xp = 1:n_quant_plt
    for xd = 1:n_quant_div
        x_mn(xp,xd) = mean(tbl.(xvar_plt)(xplt_idx==xp & xdiv_idx==xd));
        y_mn(xp,xd) = mean(tbl.(yvar)(xplt_idx==xp & xdiv_idx==xd));
        y_se(xp,xd) = std(tbl.(yvar)(xplt_idx==xp & xdiv_idx==xd))./...
            sqrt(sum(xplt_idx==xp & xdiv_idx==xd));
    end
end

%% Plot partial dependence curves
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_plt num2str(n_quant_plt) ...
    '_by_' xvar_div num2str(n_quant_div) '_quant_line_interaction'];
% figure('Name',fig_name);
hold on;
xdiv_lines = gobjects([n_quant_div 1]);
for xd = 1:n_quant_div
    xdiv_lines(xd) = errorbar(x_mn(:,xd),y_mn(:,xd),y_se(:,xd),...
        'Color',xdiv_colors(xd,:),'LineWidth',1.5,'LineStyle',xdiv_styles{xd});
end
xlabel(xplt_lab);
ylabel(yvar_lab);
title(['Interaction of ' xplt_lab{1} ' vs ' strrep(xdiv_lab{1},' Q1','')]);
set(gca,'FontSize',16);
legend(xdiv_lines,xdiv_lab,'Location','best');

end