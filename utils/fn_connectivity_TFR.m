function [conn] = fn_connectivity_TFR(tfr,metric)
%% Compute connectivity for each time-frequency point
% INPUTS:
%   tfr [ft struct] - time-frequency output of ft_freqanalysis
%       should be power for 'ampcorr', angle for 'PLV', and complex for 'coh'
%   metric [str] - name of the connectivity metric
%       {'ampcorr','PLV'} (future will add 'coh' for coherence)
% OUTPUTS:
%   conn [ft struct] - .(metric) field has [n_freq,n_time] matrix of connectivity values

if strcmp(metric,'ampcorr')
    field_name = 'powspctrm';
elseif any(strcmp(metric,{'PLV','coh'}))
    field_name = 'fourierspctrm';
else
    error(['unknown metric: ' metric]);
end

if length(tfr.label)>2; error('only ready for two channels'); end

conn = rmfield(tfr,field_name);
conn.label = [tfr.label{1} '-' tfr.label{2}];
switch metric
    case 'ampcorr'
        conn.(metric) = nan([length(tfr.freq) length(tfr.time)]);
        for f_ix = 1:length(tfr.freq)
            for t_ix = 1:length(tfr.time)
                corr = corrcoef(tfr.(field_name)(:,1,f_ix,t_ix),tfr.(field_name)(:,2,f_ix,t_ix));
                conn.(metric)(f_ix,t_ix) = corr(1,2);
            end
        end
    case 'PLV'
        % Angle diff at each time point for each trial
        angle_diff = squeeze(diff(angle(tfr.fourierspctrm),[],2));
        % Project angles to polar space, average across trials, then take magnitude of vectors
        conn.(metric) = squeeze(abs(mean(exp(1i*angle_diff),1)));
    otherwise
        error(['Unknown single trial connectivity metric: ' metric]);
end

end
