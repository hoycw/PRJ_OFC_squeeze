function [conn] = fn_connectivity_single_trial(data,metric)
%% Compute single trial connectivity
% INPUTS:
%   data [ft struct] - output of ft_preprocessing filtered data
%       should be power for 'ampcorr', angle for 'PLV', and complex for 'coh'
%   metric [str] - name of the connectivity metric
%       {'ampcorr','PLV','cohfh','cohfhl'}
%       'cohfh' = coherence using filter-Hilbert with complex data at each
%                   time point on each trial

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
    case 'cohfh'
        if ~any(imag(data.trial{1}(:))); error('coherence requires complex data!'); end
        conn = nan([n_trl 1]);
        for trl_ix = 1:n_trl
            % Fast computation of cross-spectral density and auto-spectral density at each time point:
            csd  = data.trial{trl_ix}(1,:).*conj(data.trial{trl_ix}(2,:));
            asd1 = data.trial{trl_ix}(1,:).*conj(data.trial{trl_ix}(1,:));
            asd2 = data.trial{trl_ix}(2,:).*conj(data.trial{trl_ix}(2,:));
            % Compute magnitude-squared coherence:
            %   Average across time, square the magnitude of vectors, 
            %   normalize by average power across time in 2 channels
            conn(trl_ix) = abs(mean(csd))^2 / (mean(asd1)*mean(asd2));
            
            % Per MX Cohen, this is a slower but more intuitive
            % computation of csd, which I confirmed yields differences on the
            % order of ~10^-19
%             sig1 = data.trial{trl_ix}(1,:); sig2 = data.trial{trl_ix}(2,:);
%             csd = abs(sig1).*abs(sig2).*exp(1i*(angle(sig1)-angle(sig2)));
        end
        % Originally, was computing after saving csd/asd on all trials:
%         conn = abs(mean(csd,2)).^2 ./ (mean(asd1,2).*mean(asd2,2));
    otherwise
        error(['Unknown single trial connectivity metric: ' metric]);
end

end
