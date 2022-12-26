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
    case 'logz' % log transform then z-score
        if any(pred(:)<0)
            shift_val = -min(pred) + 0.001;
            log_pred = log(pred + shift_val) - shift_val;
        else
            log_pred = log(pred);
        end
        predictor = (log_pred-nanmean(log_pred))./nanstd(log_pred);
    case 'fishz'
        predictor = atanh(pred);
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