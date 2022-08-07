%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
clc
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};
man_trl_rej_ix = {[], [71 72], [], [27 28 79 80 86 87 97 98 102 103 128 139 140 148 149 150]};
% These are true trl_ix after accounting for bad behavioral trials
%   tossed trials are based on maxabs using:
%   outlier_std_thresh = 2
%   evnt_lab = 'full' (so covering any data points included in TFR calculations
%   trl_lim = [-2 2.75];
%       2 hz min freq with 7 cycle wavelet = 1.75s added to trl_lim
%       baseline -250:50 ms means -2, then RT+1+1.75 = 2.75

plot_it  = 0;           % 0/1 whether to creat plots
evnt_lab = 'stim';      % event to time-lock segmented data {'stim' | 'dec' | 'full' = [lim(1)+S lim(2)+D]}
trl_lim = [-2 9];       % cut trial data (in sec) relative to event
new_srate = 1000;
trl_lim_samp = trl_lim.*new_srate;
n_choice_trl = 75;

rt_thresh = 6;              % hard threshold in sec for outlier RTs
amp_thresh = 100;           % hard threshold in voltage for outlier artifact trials
outlier_std_thresh = 2;

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
        
        % Loading neural data
        load(nrl_fname);
        
        % PFC01 - put it back into signal.
        if strcmp(SBJs{s},'PFC01') && ~exist('signal','var')
            signal(1,:)=datar(:,1)';        % 1 = +2 -0      % 2 = +1 -3
            signal(2,:)=datar(:,3)';        % 3 = +9 -8    % 4 = +11 -10
        end
        
        % Synchronise data- CWH: resample to 1 kHz
        % Set the resample rate of PCS data to be the same as the sampling rate in the Matlab /
        % Biopacs
        srate = SyncDetails(s,7);%422;  % Sample Rate
        n_samp = size(signal,2);
        srorig = srate;
        len_s  = n_samp./srorig;
        dt = 1./srate;
        orig_time_ax = [0:dt:len_s]; orig_time_ax(end)=[];
        nrl_resamp = [];
        for ch_ix = 1:2
            nrl_resamp(ch_ix,:) = resample(signal(ch_ix,:),orig_time_ax,new_srate);
        end
        
        %% Plot PSDs of the data
        if plot_it
            [pxx1,f] = pwelch(nrl_resamp(1,:),new_srate,0,[1:100],new_srate);
            [pxx2,f] = pwelch(nrl_resamp(2,:),new_srate,0,[1:100],new_srate);
            
            fig_name = [SBJs{s} '_run_' num2str(r_ix) '_PSDs'];
            figure; loglog(f,pxx1); hold on; loglog(f,pxx2); box off
            xticks([1 4 8 12 20 30:10:100]);
            legend('Signal 1- LFP','Signal 2 - PFC'); legend boxoff
            title(fig_name);
            
            fig_fname = [fig_dir 'PSDs/' fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
        
        %% Find and clean up the PCS data
        %     figure; plot(dtrs1(1,:));
        % PCSstr=input('Input timing of the final pulse   ');
        nrl_trim = nrl_resamp(:,PCS_start:end);
        
        load(bhv_fname);
        
        % Matlab timer baseline time is at the beginning of the first synchronisation acquisition
        matlab_bsln_time = syncR1.strtm;
        % Find and sync the Biopacs / Matlab data.
        if plot_it
            figure; plot(syncR1.data(:,3)); title([SBJs{s} ' block ' num2str(r_ix) ' Biopac data']);
        end
        % Biopstr=input('Input timing of the final pulse');
        
        %% Compile all neural and behavioral data       
        % Create Fieldtrip structure
        run_data{r_ix} = struct;
        % data.fsample=new_srate;
        run_data{r_ix}.label{1,1} = 'LFP';
        run_data{r_ix}.label{2,1} = sbj_pfc_roi{s};
        run_data{r_ix}.fsample    = new_srate;
        
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
        run_data{r_ix}.time  = cell([n_choice_trl 1]);
        run_data{r_ix}.trial = cell([n_choice_trl 1]);
%         run_data{r_ix}.sampleinfo = nan([n_choice_trl 2]);
        for trl_ix = choice_trl_ix
            % Add behavioral data (loop through to cover empty
            if ~isempty(bhv_raw(trl_ix).key)
                key{r_ix}(trl_ix) = bhv_raw(trl_ix).key;
            end
            
            % Find event timestamps
            trial_onset = bhv_raw(trl_ix).startStim-PCS_comb_bhv_start;  %1.5 seconds until motor mapping revealed.
            if bhv_raw(trl_ix).Yestrial==1
                choice_onset = bhv_raw(trl_ix).YesChoice-PCS_comb_bhv_start;
            else
                choice_onset = bhv_raw(trl_ix).NoChoice-PCS_comb_bhv_start;
            end
            rt{r_ix}(trl_ix) = choice_onset-trial_onset;
            
            % Convert to neural timestamp
            if strcmp(evnt_lab,'full')
                start_samp = round(trial_onset*new_srate) + trl_lim_samp(1);
                end_samp   = round(choice_onset*new_srate) + trl_lim_samp(2);
            else
                if strcmp(evnt_lab,'stim')
                    event_samp = round(trial_onset*new_srate);
                elseif strcmp(evnt_lab,'dec')
                    event_samp = round(choice_onset*new_srate);
                else
                    error('not ready for event besides stim or dec');
                end
                start_samp = event_samp + trl_lim_samp(1);
                end_samp   = event_samp + trl_lim_samp(2);
            end
            
            % Segment neural data
            if all(~isnan([start_samp end_samp])) % skip choice trials with no response
                %   Remember this is dependent on sample rate and should be in seconds.
                run_data{r_ix}.trial{trl_ix} = nrl_trim(:,start_samp:end_samp);
                run_data{r_ix}.time{trl_ix} = ([start_samp:end_samp]-start_samp+trl_lim_samp(1))./new_srate;
                run_data{r_ix}.sampleinfo(trl_ix,:) = [start_samp end_samp];
%                 run_data.trialinfo(trl_ix,1)=trls1(trl_ix).effortIx;
%                 run_data.trialinfo(trl_ix,2)=trls1(trl_ix).stakeIx;
            end
        end
        
        % Create previous trial regressors
        run_effort_prv{r_ix} = run_effort{r_ix};
        run_effort_prv{r_ix}(end) = [];                     % remove last trial
        run_effort_prv{r_ix} = [nan; run_effort_prv{r_ix}]; % add nan to start, shifting by 1
        
        run_stake_prv{r_ix} = run_stake{r_ix};
        run_stake_prv{r_ix}(end) = [];
        run_stake_prv{r_ix} = [nan; run_stake_prv{r_ix}];
        
        run_decision_prv{r_ix} = run_decision{r_ix};
        run_decision_prv{r_ix}(end) = [];
        run_decision_prv{r_ix} = [nan; run_decision_prv{r_ix}];
        
        clear syncR1 syncR2 result signal nrl_resamp nrl_trim orig_time_ax power
%         close all
    end
    
    %% Combine data across blocks
    % Combine behavioral data
    bhv = struct;
    bhv.key          = [key{1}; key{2}];
    bhv.rt           = [rt{1}; rt{2}];
    bhv.blk          = [run_blk_ix{1}; run_blk_ix{2}];
    bhv.blk_trl      = [run_trl_ix{1}; run_trl_ix{2}];
    bhv.run_trl      = [run_all_trl_ix{1}; run_all_trl_ix{2}];
    bhv.trl          = [run_all_trl_ix{1}; run_all_trl_ix{2}+run_all_trl_ix{1}(end)];
    bhv.run          = [ones(size(run_blk_ix{1})); ones(size(run_blk_ix{2}))*2];
    bhv.effort       = [run_effort{1}; run_effort{2}];
    bhv.stake        = [run_stake{1}; run_stake{2}];
    bhv.decision     = [run_decision{1}; run_decision{2}];
    bhv.effort_prv   = [run_effort_prv{1}; run_effort_prv{2}];
    bhv.stake_prv    = [run_stake_prv{1}; run_stake_prv{2}];
    bhv.decision_prv = [run_decision_prv{1}; run_decision_prv{2}];
    
    % Combine neural data
    data = run_data{1};
    for trl_ix=1:n_choice_trl
        data.trial{trl_ix+n_choice_trl} = run_data{2}.trial{trl_ix};
        data.time{trl_ix+n_choice_trl}  = run_data{2}.time{trl_ix};
        data.sampleinfo(trl_ix+n_choice_trl,:) = run_data{2}.sampleinfo(trl_ix,:);
    end
    
    %% Plot RTs
    if plot_it
        % now using hard threshold
        % rt_thresh = nanmean(bhv.rt) + nanstd(bhv.rt)*outlier_std_thresh;
        fig_name = [SBJs{s} '_RT_histogram'];
        figure('name',fig_name); hold on;
        histogram(bhv.rt,25);
        line([nanmean(bhv.rt) nanmean(bhv.rt)],ylim,'Color','k','LineWidth',3);
        line([rt_thresh rt_thresh],ylim,'Color','r','LineWidth',3);
        xlabel('Reaction Time (s)');
        title(SBJs{s});
        set(gca,'FontSize',16);
        
        fig_fname = [fig_dir 'RT_hist/' fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
    
    %% Remove trials with bad behavioral data
    bhv.good_key_ix = find(bhv.key==37 | bhv.key==39);
    bhv.empty_ix    = find(isnan(bhv.key));
    bhv.bad_resp_ix = find(~isnan(bhv.key) & bhv.key~=37 & bhv.key~=39);
    bhv.bad_rt_ix   = find(bhv.rt > rt_thresh);
    bad_ix = unique([bhv.empty_ix,bhv.bad_resp_ix,bhv.bad_rt_ix]);
    if ~isempty(bad_ix)
        fprintf(2,'\t%s: Removing %d empty, %d bad response, and %d bad RT trials...\n',...
            SBJs{s},length(bhv.empty_ix),length(bhv.bad_resp_ix),length(bhv.bad_rt_ix));
        % Remove from neural
        data.time(bad_ix)  = [];
        data.trial(bad_ix) = [];
        
        % Remove from behavior
        bhv_fields = fieldnames(bhv);
        for f = 1:length(bhv_fields)
            if length(bhv.(bhv_fields{f}))==n_choice_trl*2 && ~contains(bhv_fields{f},'_ix')
                bhv.(bhv_fields{f})(bad_ix) = [];
            end
        end
    end
    
    %% BEHAVIORAL MODELLING %%
    % Fit the behavior. Minimise the difference between the probability and the decision
    decisionfun=@(p) norm( (exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))) ./ ...
        (exp(p(1)) + exp(p(1)*(bhv.stake-(p(2)*(bhv.effort).^2))))) - bhv.decision);
    [par, fit]=fminsearch(decisionfun, [1,1]);
    
    SV_fn    = @(k) bhv.stake-(k*(bhv.effort).^2);
    EFF_fn   = @(k) (k*(bhv.effort).^2);
    bhv.SV   = SV_fn(par(2));
    bhv.EFFs = EFF_fn(par(2));
    bhv.p_accept = (exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2))) ./...
        (exp(par(1)) + exp(par(1)*(bhv.stake-(par(2)*(bhv.effort).^2)))));
    
    % Creat previous trial predictors (respecting tossed trials)
    bhv.SV_prv       = nan(size(bhv.trl));
    bhv.EFFs_prv     = nan(size(bhv.trl));
    bhv.p_accept_prv = nan(size(bhv.trl));
    firsttrl_ix = find(bhv.run_trl==1);
    for t_ix = 1:length(bhv.trl)
        prv_ix = find(bhv.trl==bhv.trl(t_ix)-1);
        if ~isempty(prv_ix) && ~any(t_ix==firsttrl_ix)  % only for previous trial and not first trial of the run
            bhv.SV_prv(t_ix)       = bhv.SV(prv_ix);
            bhv.EFFs_prv(t_ix)     = bhv.EFFs(prv_ix);
            bhv.p_accept_prv(t_ix) = bhv.p_accept(prv_ix);
        end
    end
    
    if plot_it
        figure;
        subplot(3,1,1);
        scatter(bhv.stake,bhv.SV)
        xlabel('Apples'); ylabel('SV');
        subplot(3,1,2);
        scatter(bhv.effort,bhv.SV)
        xlabel('Effort - Proportion max'); ylabel('SV');
        subplot(3,1,3);
        scatter(bhv.SV,bhv.p_accept)
        xlabel('SV'); ylabel('Probabilty of acceptance');
        
        % Plot a regression line %
%         SVtmp = bhv.SV;
%         
%         figure; scatter(SVtmp,ProbAccept,'.')
%         
%         sigparam=sigm_fit(SVtmp,ProbAccept,1)
%         
%         [p,S,mu]  = polyfit(SVtmp,ProbAccept,12)
%         PotValues=[-10:0.1:10];
%         Pout=polyval(p,PotValues);
%         
%         figure; plot(PotValues,Pout)
        
        % Test behavioral modelling.
        figure;
        subplot(2,1,1);
        scatter(bhv.stake, bhv.SV);
        title('Reward')
        subplot(2,1,2);
        scatter(bhv.effort, bhv.SV);
        title('Effort');
    end
    
    % Regression
    % Create regression
%     RG    = [ones(size(bhv.stake)),bhv.stake,bhv.effort];
%     betas = regress(bhv.SV,RG); %% FIELDTRIP TF ANALYSIS
    
    % Start by ignoring the autonomic data and just do a correlation of the
    % neural features by the decision features.
    
    %% Data Visualization for Trial Rejection
    if plot_it
        data_plot = rmfield(data,'sampleinfo');
%         cfgpp = [];
%         cfgpp.lpfilter = 'yes';
%         cfgpp.lpfreq   = 50;
%         cfgpp.hpfilter = 'yes';
%         cfgpp.hpfreq   = 0.5;
%         cfgpp.hpfiltord = 4;
%         data_pp = ft_preprocessing(cfgpp,data_plot);
        
        cfgp = []; cfgp.viewmode = 'vertical';
        db_out = ft_databrowser(cfgp, data_plot);
        
        cfgr = []; cfgr.method = 'summary';
        reject_out = ft_rejectvisual(cfgr,data_plot);
    end
    
    %% Find neural outliers
    % Get standard deviation and normalize per channel
    trial_abs = nan([size(data.trial,1) length(data.label)]);
    for trl_ix = 1:size(data.trial,1)
        tmp = data.trial{trl_ix};
        for ch_ix = 1:2
            trial_abs(trl_ix,ch_ix)  = max(abs(tmp(ch_ix,:)));
        end
%         trialmean(trl_ix,:) = mean(abs(tmp),2);    % why are you taking abs here?
    end
    trial_abs_norm = trial_abs./nanmean(trial_abs);
    
    % figure;
    % scatter(1:length(trial_abs_norm),squeeze(trial_abs_norm(:,1)));
    
    % Find trials with outlier standard deviations
    var_reject_ix = [];
    for ch_ix = 1:2
        var_reject_ix = [var_reject_ix; find(trial_abs_norm(:,ch_ix)>outlier_std_thresh)];
    end
    var_reject_ix = unique(var_reject_ix);
    
    %% Reject neural outliers
    % Rejection based on criteria listed in parameter definition
    outlier_ix = nan(size(man_trl_rej_ix{s}));
    for t = 1:numel(man_trl_rej_ix{s})
        outlier_ix(t) = find(bhv.trl==man_trl_rej_ix{s}(t));
    end
    if isempty(man_trl_rej_ix{s})
        fprintf('%s: No bad trials identified.\n',SBJs{s});
    else
        fprintf(2,'%s: Removing %d trials!\n',SBJs{s},length(man_trl_rej_ix{s}));
    end
    
    % Remove outlier trials
    data.trial(outlier_ix) = [];
    data.time(outlier_ix)  = [];
%     data.trialbl(bad_trl_ix) = [];
%     data.sampleinfo(bad_trl_ix,:) = [];
    
    bhv.reject_ix = outlier_ix;
    bhv_fields = fieldnames(bhv);
    n_good_bhv = length(bhv.trl);
    for f = 1:length(bhv_fields)
        if length(bhv.(bhv_fields{f}))==n_good_bhv && ~contains(bhv_fields{f},'_ix')
            bhv.(bhv_fields{f})(outlier_ix) = [];
        end
    end
    
    %% SAVE data
    preproc_name = [FileDetails{s,1} '_' evnt_lab]; % Patient code _ event
    sbj_data = struct;
    sbj_data.desp    = preproc_name;
    sbj_data.ts      = data;
%     sbj_data.TF      = freqout;
%     sbj_data.TFbl    = freqoutbl;
    sbj_data.bhv     = bhv;
    sbj_data.mdl.par = par;
    sbj_data.mdl.fit = fit;
    
    out_dir = [prj_dir 'data/' SBJs{s} '/'];
    out_fname = [out_dir preproc_name '_preproc.mat'];
    fprintf('Saving %s\n',out_fname);
    save(out_fname,'sbj_data');
    
end