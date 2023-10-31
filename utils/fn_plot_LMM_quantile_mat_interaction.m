function fn_plot_LMM_quantile_mat_interaction(tbl,xvar_plt,xvar_div,yvar,n_quantiles,varargin)
% Matrix plots of yvar as a function of xvar_plt, discretized by xvar_div
% (one line per value of xvar_div)

% if ~contains(xvar_prv,'prv') || ~contains(xvar_cur,'cur')
%     error('input xvar must have cur/prv for gratton violin plot');
% end
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'continuous')
            continuous_data = varargin{v+1};
        elseif strcmp(varargin{v},'logistic_fit')
            logistic_fit = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end
if ~exist('continuous_data','var'); continuous_data = 0; end
if ~exist('logistic_fit','var');   logistic_fit = 0; end

%% Set up and check vaeriables
scat_sz  = 50;
x_fudge = 0.2;

sbj_ns = unique(tbl.sbj_n);
if ~all(sbj_ns'==[1 2 3 4]); error('SBJs in tbl mismatch'); end
[xplt_lab, ~] = fn_get_label_styles(xvar_plt,1);
[xdiv_lab, ~] = fn_get_label_styles(xvar_div,1);
[yvar_lab, ~] = fn_get_label_styles(yvar,1);

%% Discretize predictors into quantiles
xdiv_idx = zeros(size(tbl.sbj_n));
xplt_idx = zeros(size(tbl.sbj_n));
for s = 1:length(sbj_ns)
    % Get bin edges
    if continuous_data
        div_qs = quantile(tbl.(xvar_div)(tbl.sbj_n==s),n_quantiles);
        plt_qs = quantile(tbl.(xvar_div)(tbl.sbj_n==s),n_quantiles);
    else
        div_qs = unique(tbl.(xvar_div)(tbl.sbj_n==s));
        plt_qs = unique(tbl.(xvar_div)(tbl.sbj_n==s));
        if length(div_qs)~=n_quantiles
            warning(['requested ' num2str(n_quantiles) ' quantiles but ' num2str(length(div_qs)) ' are in the data']);
        end
        if length(plt_qs)~=n_quantiles
            warning(['requested ' num2str(n_quantiles) ' quantiles but ' num2str(length(plt_qs)) ' are in the data']);
        end
    end
    div_edges = nan([n_quantiles-1 1]);
    plt_edges = nan([n_quantiles-1 1]);
    for i = 1:n_quantiles-1
        div_edges(i) = mean(div_qs(i:i+1));
        plt_edges(i) = mean(plt_qs(i:i+1));
    end
    
    % Assign to quantiles
    for i = 1:n_quantiles
        if i==1                 % first quantile
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<div_edges(i)) = i;
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)<plt_edges(i)) = i;
        elseif i==n_quantiles   % last quantile
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)>=div_edges(i-1)) = i;
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)>=plt_edges(i-1)) = i;
        else                    % middle quantiles
            xdiv_idx(tbl.sbj_n==s & tbl.(xvar_div)<div_edges(i) &...
                tbl.(xvar_div)>=div_edges(i-1)) = i;
            xplt_idx(tbl.sbj_n==s & tbl.(xvar_plt)<plt_edges(i) &...
                tbl.(xvar_plt)>=plt_edges(i-1)) = i;
        end
    end
end

%% Average within bins
y_mn = nan([n_quantiles n_quantiles]);
for xp = 1:n_quantiles
    for xd = 1:n_quantiles
        y_mn(xp,xd) = mean(tbl.(yvar)(xdiv_idx==xp & xplt_idx==xd));
    end
end

%% Plot partial dependence curves
fig_name = ['GRP_TFR_LMM_results_' yvar '_' xvar_plt '_by_' xvar_div '_quant_matrix_interaction'];
figure('Name',fig_name); hold on;
imagesc(1:n_quantiles,1:n_quantiles,y_mn);
xlabel(xdiv_lab);
xlim([0.5 n_quantiles+0.5]);
xticks(1:n_quantiles);
xticklabels(num2str(div_qs,'%.2f'));
ylabel(xplt_lab);
ylim([0.5 n_quantiles+0.5]);
yticks(1:n_quantiles);
yticklabels(num2str(plt_qs,'%.2f'));
set(gca,'YDir','normal');
cbar = colorbar;
cbar.Label.String = yvar_lab;
title(['Interaction of ' xplt_lab{1} ' vs ' xdiv_lab{1}]);
set(gca,'FontSize',16);

end