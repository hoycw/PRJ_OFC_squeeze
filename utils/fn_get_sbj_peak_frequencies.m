function [theta_peak_freqs, betalo_peak_freqs, betahi_peak_freqs] = fn_get_sbj_peak_frequencies(SBJs,an_id)
%% Returns the center (peak) frequency for theta and low beta for PFC squeeze SBJs
%   based on the analysis parameters; if no SBJ specific peak, then -1
%   assigned to indicate a default setting should be used
% INPUTS:
%   SBJs [cell array] - strings of SBJ names
%       should be {'PFC03','PFC04','PFC05','PFC01'}
%   an_id [str] - ID string for TFR analysis
% OUTPUTS:
%   theta_peak_freqs [2x4 float array] - theta peak frequency for each SBJ and channel
%       row 1 is basal ganglia, row 2 is PFC
%   betalo_peak_freqs [2x4 float array] - low beta peak frequency for each SBJ and channel
%       row 1 is basal ganglia, row 2 is PFC

if ~all(strcmp(SBJs,{'PFC03','PFC04','PFC05','PFC01'}))
    error('wrong SBJ order!');
end

betahi_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
if contains(an_id,'_S')
    if any(strcmp(an_id,{'TFRmth_S1t2_zS1t0_f2t40','simon'}))  % equivalent to Simon
        theta_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
        betalo_peak_freqs  = [-1 -1 -1 -1; 10 17 13 12];
    elseif any(strcmp(an_id,{'TFRmth_S1t2_zS8t0_f2t40','TFRmth_S1t2_madS8t0_f2t40'}))
        theta_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
        betalo_peak_freqs  = [11 13 17 13; 12 16 12 12]; % D-locked: [-1 -1 -1 -1; 11 19 12 12];
    elseif strcmp(an_id,'TFRw_S25t2_dbS25t05_fl2t40_c7')
        theta_peak_freqs = [2.5 3.5 3.5 6; 5 4 3.5 3];
        betalo_peak_freqs  = [-1 -1 -1 -1; -1 -1 -1 -1];
    else
        theta_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
        betalo_peak_freqs = [-1 -1 -1 -1; 10 17 13 12];
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
%         theta_peak_freqs = [2.5 3 3.5 5; 5 3.5 3.5 3]; % BG (row1) then PFC (row2)
%         betalo_peak_freqs  = [17 22 14 14; 13 17 11 15];
        %     betahi_cf  = % no SBJ-specific peaks;
        %     betalo_peak_freqs = ones([2 numel(SBJs)])*-1;
    end
elseif contains(an_id,'_D')
    if strcmp(an_id,'TFRw_D1t1_dbS25t05_fl2t40_c7')
        theta_peak_freqs = [3.5 3 4.5 6; 3 3 4.5 3];
        betalo_peak_freqs  = [13 20 17 18; -1 -1 -1 -1];
    elseif strcmp(an_id,'TFRmth_D1t1_zS1t0_f2t40')  % equivalent to Simon
        theta_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
        betalo_peak_freqs  = [-1 -1 -1 -1; -1 -1 -1 -1];%10 17 13 12];
    elseif any(strcmp(an_id,{'TFRmth_D1t1_zS8t0_f2t40','TFRmth_D1t1_madS8t0_f2t40'}))
        theta_peak_freqs = [-1 -1 -1 -1; -1 -1 -1 -1];
        betalo_peak_freqs  = [-1 -1 -1 -1; 11 19 12 12];
    else
        % As of GRP_TFR_peak_find on 0.25-1.5 from 7/5/22
        %   Use -1 for canonical bands
        theta_peak_freqs = [3.5 3 3.5 7; 4.5 3 3.5 3];
        %     betalo_peak_freqs  = [17 22 14 14; 13 17 11 15];
        betalo_peak_freqs = [20 17 11 24; -1 -1 -1 -1];
    end
end

end