function [lme_full,pred_pval] = fn_run_LMM_full_permutation_models(tbl,predictors,full_formula,fit_method,nboots)
%% Run LMM with null permutation testing for specific predictors
% Use two-sided test

% Run true model and significance testing
tic;
pred_str = [repmat('%s, ',1,length(predictors)-1) '%s'];
fprintf(['Running LMM permutation testing:\n\tModel: %s\n\tPredictors: ' pred_str '\n'],...
    full_formula,predictors{:});
lme_full = fitlme(tbl,full_formula,'FitMethod',fit_method);

% pred_lrt   = cell(size(predictors));
% for p = 1:length(predictors)
%     null_formula = strrep(full_formula,[predictors{p} ' + '],'');
%     pred_lme = fitlme(tbl,null_formula,'FitMethod',fit_method);
%     pred_lrt{p} = compare(pred_lme,lme_full,'CheckNesting',true);
% end

% Run permutation model
null_pred_coefs = nan([length(predictors) nboots]);
parfor b_ix = 1:nboots
    null_pred_coefs(:,b_ix) = fn_run_LMM_null_permutation_iteration(...
        tbl,predictors,full_formula,fit_method);
end
fprintf('\nFinished null permutations, time elapsed = %.2f min\n\n',toc/60);

% Compute and report p values from null distribution of coefficients
fprintf('LMM Permutation results:\n');
fprintf('\tModel: %s\n',full_formula);
pred_pval = nan(size(predictors));
for p = 1:length(predictors)
    coef_ix = strcmp(lme_full.CoefficientNames,predictors{p});
    pred_pval(p) = sum(abs(null_pred_coefs(p,:))>=abs(lme_full.Coefficients.Estimate(coef_ix)))/nboots;
    fprintf('\t%s: coefficient = %.4f, p = %.4f\n',predictors{p},...
        lme_full.Coefficients.Estimate(coef_ix),pred_pval(p));
end

end