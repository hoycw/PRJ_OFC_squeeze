function [null_coefs, null_pvals] = fn_run_LMM_null_permutation(tbl,predictors,lmm_formula,fit_method)
% INPUTS:

% Permute trials separately for each subject and run, but
%   keep same trl_n permutation mapping across hemispheres
null_idx = nan(size(tbl.(predictors{1})));
null_sbjs = unique(tbl.sbj_n);
for s = 1:length(null_sbjs)
    sbj_pred = tbl.(predictors{1})(tbl.sbj_n==null_sbjs(s));
    null_idx(tbl.sbj_n==null_sbjs(s)) = randperm(length(sbj_pred));
end
if any(isnan(null_idx)); error('NaNs in null predictor!'); end

% Run null models with same permuted index for one predictor at a time
pred_lme  = cell(size(predictors));
null_coefs = nan(size(predictors));
null_pvals = nan(size(predictors));
for p = 1:length(predictors)
    % Permute predictor of interest only
    null_tbl = tbl;
    null_tbl.(predictors{p}) = tbl.(predictors{p})(null_idx);
    
    % Run full null model
    null_mdl = fitlme(null_tbl,lmm_formula,'FitMethod',fit_method);
    coef_ix = strcmp(null_mdl.CoefficientNames,predictors{p});
    null_coefs(p) = null_mdl.Coefficients.Estimate(coef_ix);
    
    % Run reduced null model
    null_bsln_formula = strrep(lmm_formula,[predictors{p} ' + '],'');
    pred_lme{p}  = fitlme(null_tbl,null_bsln_formula,'FitMethod',fit_method);
    pred_lrt = compare(pred_lme{p},null_mdl,'CheckNesting',true);
    null_pvals(p) = pred_lrt.pValue(2);
end

end