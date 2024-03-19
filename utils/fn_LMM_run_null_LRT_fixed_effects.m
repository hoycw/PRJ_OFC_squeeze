function [null_mdls, comparisons] = fn_LMM_run_null_LRT_fixed_effects(tbl,full_mdl,preds,fit_method)
%% Run null models removing each fixed effet from full LMM and return stats from likelihood ratio tests

null_mdls   = cell(size(preds));
comparisons = cell(size(preds));
for p = 1:length(preds)
    null_formula = [char(full_mdl.Formula) '- ' preds{p}];
    null_mdls{p} = fitlme(tbl,null_formula,'FitMethod',fit_method);
    comparisons{p} = compare(null_mdls{p},full_mdl,'CheckNesting',true);
end

end