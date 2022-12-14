%% Check raw behavioral data
bhv_name = 'Apples_001_01PFC05_1219Bl1ABB.mat';% 'Apples_001_01PFC05_1219Bl3ABA.mat';%'Apples_001_01PFC05_1219Bl2ABA.mat';%

rt_thresh = 6;              % hard threshold in sec for outlier RTs
outlier_std_thresh = 2;
n_choice_trl = 75;

%% Load data
data_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/data/PFC05/Session2_withStimulation/';
bhv_fname = [data_dir bhv_name];
load(bhv_fname);
bhv_raw = result.data;

%% QA checks on raw behavioral data
choice_trl_ix   = find(~isnan([bhv_raw.yeslocation]));
if numel(choice_trl_ix)~=n_choice_trl || ~all(choice_trl_ix==1:n_choice_trl)
    error('Why are the first 75 trials not choice trials?\n');
end
if length(bhv_raw)~=n_choice_trl+10
    warning('\t%s has different number of trials! (%d)\n',bhv_name,length(bhv_raw));
end

blk_n = unique([bhv_raw.block]);
for b = 1:length(blk_n)
    fprintf('block %d length = %d\n',blk_n(b),sum([bhv_raw.block]==blk_n(b)));
end

% Segment neural data based on event timestamps
run_effort     = [bhv_raw(choice_trl_ix).effort]';  %.effortIx
bhv.key = nan(size(run_effort));
bhv.rt  = nan(size(run_effort));
for trl_ix = choice_trl_ix
    % Add behavioral data (loop through to cover empty
    if ~isempty(bhv_raw(trl_ix).key)
        bhv.key(trl_ix) = bhv_raw(trl_ix).key;
    end
    
    % Find event timestamps
    trial_onset = bhv_raw(trl_ix).startStim;  %1.5 seconds until motor mapping revealed.
    if bhv_raw(trl_ix).Yestrial==1
        choice_onset = bhv_raw(trl_ix).YesChoice;
    else
        choice_onset = bhv_raw(trl_ix).NoChoice;
    end
    bhv.rt(trl_ix) = choice_onset-trial_onset;
    
end

bhv.good_key_ix = find(bhv.key==37 | bhv.key==39);
bhv.empty_ix    = find(isnan(bhv.key));
bhv.bad_resp_ix = find(~isnan(bhv.key) & bhv.key~=37 & bhv.key~=39);
bhv.bad_rt_ix   = find(bhv.rt > rt_thresh);
bad_ix = unique([bhv.empty_ix,bhv.bad_resp_ix,bhv.bad_rt_ix]);
if ~isempty(bad_ix)
    fprintf(2,'\t%s: Removing %d empty, %d bad response, and %d bad RT trials...\n',...
        bhv_name,length(bhv.empty_ix),length(bhv.bad_resp_ix),length(bhv.bad_rt_ix));
else
    fprintf('\t%s: No bad trials, all good!\n',bhv_name);
end