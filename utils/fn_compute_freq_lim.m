function [pow_lim] = fn_compute_freq_lim(SBJs,peak_freqs,band_id)

if ~all(strcmp(SBJs,{'PFC03','PFC04','PFC05','PFC01'}))
    warning('not full list or wrong SBJ order!');
end

%% Set defaults
switch band_id
    case 'theta'
        cannonical = [4 7];
        bandwidth  = 3;
    case 'betalo'
        cannonical = [12 20];
        bandwidth  = 4;
    case 'betahi'
        cannonical = [21 35];
        bandwidth  = 4;
    case 'beta'
        cannonical = [12 35];
        bandwidth  = 3;
    otherwise
        error(['frequency band not recognized: ' band_id]);
end
% Simon originals:
% alternatives: (1)=[17,17,13,12]; (2)=[17,22,13,13];

%% Compute band limits
pow_lim  = nan(length(SBJs),2,2);
for s = 1:length(SBJs)
    for ch_ix = 1:2
        if peak_freqs(ch_ix,s)==-1
            pow_lim(s,ch_ix,:) = cannonical;
        else
            pow_lim(s,ch_ix,:) = [peak_freqs(ch_ix,s)-bandwidth/2 peak_freqs(ch_ix,s)+bandwidth/2];
        end
        % Print final parameters
%         if ch_ix==1; ch_lab = sbj_bg_roi{s}; else; ch_lab = sbj_pfc_roi{s}; end
%         fprintf('%s %s beta lo CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betalo_cf(ch_ix,s),betalo_lim(s,ch_ix,1),betalo_lim(s,ch_ix,2));
%         fprintf('%s %s beta hi CF = %.02f; BW = %.02f to %.02f\n',SBJs{s},ch_lab,betahi_cf(ch_ix,s),betahi_lim(s,ch_ix,1),betahi_lim(s,ch_ix,2));
    end
end


end