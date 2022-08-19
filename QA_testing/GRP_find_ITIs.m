%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

n_choice_trl = 75;

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
fig_dir   = [prj_dir 'results/preproc/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
fig_ftype = 'png';

%% Load parameters %%
% Columns in SqueezeSubjectSyncSummary.xlsx:
%   Patient, Block1 Neural Data, Block1 Behavioral Data, Blck1 Sync1, Blck1
%   Sync 2, Block2 Neural Data, Block1 Behavioral Data (likely typo?),
%   Blck1 Sync 1, Blck1 Sync 2, SampleRate
DataStorage  = [prj_dir 'box/Analysis/Group/Preprocess/'];
DataStorage2 = [prj_dir 'box/Analysis/Group/Preprocess/Output_files_shifted_behavior/'];

[numbers, strings, raw] = xlsread(strcat(DataStorage,'SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);
SyncDetails = numbers;

%% Process data
dec_to_trl = cell(size(SBJs));
for s=1:4
%     clearvars -except s SBJs sbj_pfc_roi evnt_lab trl_lim new_srate trl_lim_samp n_choice_trl Flnum SyncDetails FileDetails
%     close all;
    
    for r_ix = 1:2
        % Set SBJ-specific filenames and sync details
        if r_ix==1
            nrl_fname = fullfile(DataStorage,'DataFiles',FileDetails{s,2});    % Block1 Neural Data
            bhv_fname = fullfile(DataStorage,'DataFiles',FileDetails{s,3});    % Block1 Behavioral Data
            PCS_start    = SyncDetails(s,1);  % Blck1 Sync 1
            biopac_start = SyncDetails(s,2);  % Blck1 Sync 2
        elseif r_ix==2
            nrl_fname = fullfile(DataStorage,'DataFiles',FileDetails{s,6});    % Block2 Neural Data
            bhv_fname = fullfile(DataStorage,'DataFiles',FileDetails{s,7}); % Block2? behavioral data
            PCS_start    = SyncDetails(s,5);    % Blck2? Sync 1
            biopac_start = SyncDetails(s,6);  % Blck2? Sync 2
        else
            error('only 2 blocks!');
        end
        
        %% Find and clean up the PCS data
        load(bhv_fname);
        
        % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
        matlab_bsln_time = syncR1.strtm;
        
        %% Compile all neural and behavioral data       
        ad1 = biopac_start./500;
        % Have checked the timings of this and it works now.
        % !!! CWH: unclear what these timestamps are, but Simon says above it's okay...
        PCS_comb_bhv_start = matlab_bsln_time+ad1;
        
        % 
        bhv_raw = result.data;
        choice_trl_ix   = find(~isnan([bhv_raw.yeslocation]));
        if numel(choice_trl_ix)~=n_choice_trl || ~all(choice_trl_ix==1:n_choice_trl)
            error('Why are the first 75 trials not choice trials?\n');
        end
        if length(bhv_raw)~=n_choice_trl+10
            warning('\t%s run #%d has different number of trials! (%d)\n',SBJs{s},r_ix,length(bhv_raw));
        end
        
        % Compile behavioral data
%             effortlevels = [0.16 0.32 0.48 0.64 0.80];
%             stakelevels  = [1 4 7 10 13];
        run_effort{r_ix}     = [bhv_raw(choice_trl_ix).effort]';  %.effortIx
        run_stake{r_ix}      = [bhv_raw(choice_trl_ix).stake]';   %.stakeIx
        run_decision{r_ix}   = [bhv_raw(choice_trl_ix).Yestrial]';
        run_blk_ix{r_ix}     = [bhv_raw(choice_trl_ix).block]';
        run_trl_ix{r_ix}     = [bhv_raw(choice_trl_ix).trialIndex]';
        run_all_trl_ix{r_ix} = [bhv_raw(choice_trl_ix).allTrialIndex]';
        
        % Segment neural data based on event timestamps
        key{r_ix} = nan(size(run_effort{r_ix}));
        rt{r_ix}  = nan(size(run_effort{r_ix}));
        trial_onsets{r_ix}   = nan(size(run_effort{r_ix}));
        choice_onsets{r_ix}  = nan(size(run_effort{r_ix}));
        run_data{r_ix}.time  = cell([n_choice_trl 1]);
        run_data{r_ix}.trial = cell([n_choice_trl 1]);
%         run_data{r_ix}.sampleinfo = nan([n_choice_trl 2]);
        for trl_ix = choice_trl_ix
            % Add behavioral data (loop through to cover empty
            if ~isempty(bhv_raw(trl_ix).key)
                key{r_ix}(trl_ix) = bhv_raw(trl_ix).key;
            end
            
            % Find event timestamps
            trial_onsets{r_ix}(trl_ix) = bhv_raw(trl_ix).startStim-PCS_comb_bhv_start;  %1.5 seconds until motor mapping revealed.
            if bhv_raw(trl_ix).Yestrial==1
                choice_onsets{r_ix}(trl_ix) = bhv_raw(trl_ix).YesChoice-PCS_comb_bhv_start;
            else
                choice_onsets{r_ix}(trl_ix) = bhv_raw(trl_ix).NoChoice-PCS_comb_bhv_start;
            end
%             rt{r_ix}(trl_ix) = choice_onset-trial_onset;
            
        end
        clear syncR1 syncR2 result signal nrl_resamp nrl_trim orig_time_ax power
%         close all
    end
    
    %% Compute ITIs
    trl_on = [trial_onsets{1}; trial_onsets{2}];
    dec_on = [choice_onsets{1}; choice_onsets{2}];
    dec_on_prv = [nan; dec_on(1:end-1)];
    dec_to_trl{s} = trl_on-dec_on_prv;
    
    blk_trl      = [run_trl_ix{1}; run_trl_ix{2}];
    dec_to_trl{s}(blk_trl==1) = [];
end

%% Plot ITI distributions
figure;
for s = 1:4
    subplot(2,2,s);
    histogram(dec_to_trl{s},25);
    xlabel('Time from decision to trial onset (s)');
    title(SBJs{s});
    set(gca,'FontSize',16);
end

fig_fname = [fig_dir 'SBJ_ITI_hist.' fig_ftype];
fprintf('Saving %s\n',fig_fname);
saveas(gcf,fig_fname);
