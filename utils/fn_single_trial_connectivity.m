function [conn] = fn_single_trial_connectivity(data,metric)
%% Compute single trial connectivity
% INPUTS:
%   data [ft struct] - output of ft_preprocessing filtered data
%       should be power for 'ampcorr', angle for 'PLV', and complex for 'coh'
%   metric [str] - name of the connectivity metric
%       {'ampcorr','PLV'} (future will add 'coh' for coherence)

n_trl = length(data.trial);
if length(data.label)>2; error('only ready for two channels'); end

switch metric
    case 'ampcorr'
        conn = nan([n_trl 1]);
        for trl_ix = 1:n_trl
            corr = corrcoef(data.trial{trl_ix}(1,:),data.trial{trl_ix}(2,:));
            conn(trl_ix) = corr(1,2);
        end
    case 'PLV'
        % Angle diff at each time point for each trial
        angle_diff = nan([n_trl length(data.time{1})]);
        for trl_ix = 1:n_trl
            angle_diff(trl_ix,:) = data.trial{trl_ix}(1,:)-data.trial{trl_ix}(2,:);
        end
        % Project angles to polar space, average across time, then take magnitude of vectors
        conn = abs(mean(exp(1i*angle_diff),2));
    otherwise
        error(['Unknown single trial connectivity metric: ' metric]);
end

end