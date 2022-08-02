function [predictor] = fn_normalize_predictor(pred,norm_method)
%% Normalize a model predictor (or not)
% INPUTS:
%   pred [vector] - predictor to be normalized
%   norm_method [str] - method to use to normalize vector
% OUTPUTS:
%   predictor [vector] - normalized (or not) predictor

if size(pred,2)>1 || numel(size(pred))>2
    error('pred is not a column vector!');
end

switch norm_method
    case 'zscore'
        predictor = (pred-nanmean(pred))./nanstd(pred);
    case 'demean'
        predictor = pred-nanmean(pred);
    case 'minmax'
        predictor = (pred-min(pred))./(max(pred)-min(pred));
    case 'none'
        predictor = pred;
    otherwise
        error(['Unknown normalization method: ' norm_method]);
end

end