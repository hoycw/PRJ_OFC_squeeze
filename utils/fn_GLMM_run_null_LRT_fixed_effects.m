function [null_mdls, comparisons] = fn_GLMM_run_null_LRT_fixed_effects(tbl,full_mdl,preds,fit_method)
%% Run null models removing each fixed effet from full LMM and return stats from likelihood ratio tests

null_mdls   = cell(size(preds));
comparisons = cell(size(preds));
for p = 1:length(preds)
    formula = char(full_mdl.Formula);
    if contains(formula,[' + ' preds{p} ' '])
        % Remove that term directly
        null_formula = strrep(formula,['+ ' preds{p}],'');
    else
        % Try to remove that term just before the first random effects term
        re_pos = strfind(formula,'+ (');
        if ~isempty(re_pos)
            null_formula = insertBefore(formula,re_pos(1),['- ' preds{p} ' ']);
        else % no random effects, just add to end
            null_formula = [formula ' - ' preds{p}];
        end
    end
    null_mdls{p} = fitglme(tbl,null_formula,'Distribution','binomial','FitMethod',fit_method);
    comparisons{p} = compare(null_mdls{p},full_mdl,'CheckNesting',true);
end

end