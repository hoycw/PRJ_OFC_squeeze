function fn_plot_LMM_quantile_lines(SBJs,tbl,xvar,yvar,lmm,pval,n_quantiles)

%% Set up and check vaeriables
if ~all(strcmp(SBJs,{'PFC03','PFC04','PFC05','PFC01'})); error('SBJs wrong'); end
if ~all(unique(tbl.sbj_n)'==[1 2 3 4]); error('SBJs in tbl mismatch'); end

x_fudge = 0.2;
sbj_colors = [27, 158, 119;         % teal
              217, 95, 2;           % burnt orange
              117, 112, 179;        % purple
              231, 41, 138]./256;   % magenta
coef_ix = strcmp(lmm.CoefficientNames,xvar);

%% Discretize predictor into quantiles
for s = 1:length(SBJs)
    % Get bin edges
    qs = quantile(tbl.(xvar)(tbl.sbj_n==s),n_quantiles);
    edges = nan([n_quantiles-1 1]);
    for i = 1:n_quantiles-1
        edges(i) = mean(qs(i:i+1));
    end
    
    % Average within bins
    x_mn{s} = nan([n_quantiles 1]);
    y_mn{s}  = nan([n_quantiles 1]);
    y_se{s}  = nan([n_quantiles 1]);
    for i = 1:n_quantiles
        if i==1                 % first quantile
            trl_idx = tbl.sbj_n==s & tbl.(xvar)<edges(i);
        elseif i==n_quantiles   % last quantile
            trl_idx = tbl.sbj_n==s & tbl.(xvar)>=edges(i-1);
        else                    % middle quantiles
            trl_idx = tbl.sbj_n==s & tbl.(xvar)<edges(i) &...
                tbl.(xvar)>=edges(i-1);
        end
        x_mn{s}(i) = mean(tbl.(xvar)(trl_idx));
        y_mn{s}(i) = mean(tbl.(yvar)(trl_idx));
        y_se{s}(i) = std(tbl.(yvar)(trl_idx))./sqrt(sum(trl_idx));
    end
end

%% Create plot
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar '_quant_line'];
figure('Name',fig_name); hold on;
xvals = min(tbl.(xvar))-x_fudge:0.01:max(tbl.(xvar))+x_fudge;
for s = 1:length(SBJs)
    errorbar(x_mn{s},y_mn{s},y_se{s},...
        'Color',sbj_colors(s,:),'LineWidth',1.5);
end
yvals = lmm.Coefficients.Estimate(1) + xvals*lmm.Coefficients.Estimate(coef_ix);
line(xvals,yvals,'Color','k','LineWidth',4);
xlabel(xvar);
ylabel(yvar);
title(['LMM coef. = ' num2str(lmm.Coefficients.Estimate(coef_ix)) '; p = ' num2str(pval,'%.03f')]);
legend([SBJs 'Group'],'Location','best');
set(gca,'FontSize',16);


end