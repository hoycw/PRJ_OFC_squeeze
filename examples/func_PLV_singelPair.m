function PLV_across_trials = func_PLV_singlePair(hilbert1, hilbert2, starts, ends, foldername, elec1, elec2)
%Taking 2 electrodes and calculating average phase coherence over all trials

% Cutoff signals by time range
for k = 1:length(starts)
    a = starts(k);
    b = ends(k);
    ch1(k,:) = hilbert1(a:b);
    ch2(k,:) = hilbert2(a:b);
end

% Angle diff at each time point for each trial
angle_diff = ch1 - ch2;

%PLV across time/trial (one PLV for each trial, then avg/one PLV for each time point, then avg)
PLV_vector = abs(mean(exp(1i*angle_diff),2));
PLV_across_trials = mean(PLV_vector);
PLV_across_time = mean(abs(mean(exp(1i*angle_diff)))); % shouldn't we skip the outermost "mean" to get a time series?

%Save PLV across trials
cd(foldername);
save(strcat('PLVe',num2str(elec1,'%03d'),'_e',num2str(elec2,'%03d')), 'PLV_vector');
% save(strcat('angleDiff_trialXtime_e',num2str(elec1,'%03d'),'_e',num2str(elec2,'%03d')), 'angle_diff');

end




