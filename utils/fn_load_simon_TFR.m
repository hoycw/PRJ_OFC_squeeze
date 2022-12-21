function [tfr, bhv] = fn_load_simon_TFR(SBJs,toss_same_trials,man_trl_rej_ix)
%% Load original Simon TFR files
%   Option to remove bad trials

if ~all(strcmp(SBJs,{'PFC03','PFC04','PFC05','PFC01'}))
    error('wrong SBj order!');
end

%% Load Simon file details
DataStorage='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess';
DataStorage2='/Users/colinhoy/Code/PRJ_OFC_squeeze/box/Analysis/Group/Preprocess/Output_files_shifted_behavior';

[numbers, strings, raw] = xlsread(strcat(DataStorage,'/','SqueezeSubjectSyncSummary.xlsx'));
% SBJs: PFC03, PFC04, PFC05, PFC01, PMC10
if ~all(strcmp(strings(2:6,1),[SBJs {'PMC10'}]')); error('SBJ order off!'); end

FileDetails = strings(2:end,:);

%% Load data files
proc_fname = strcat(DataStorage2,'/',FileDetails{s,1},'Stimulus_Locked.mat');
fprintf('Loading %s\n',proc_fname);
load(proc_fname);

tfr = AllData.TFbl;
bhv = AllData.exp;

%% Toss outlier trials
if toss_same_trials
    load([sbj_dir SBJs{s} '_stim_preproc.mat'],'sbj_data');
    % Combine bad behavioral and neural trials
    %   (empty and bad key are already tossed in Simon data)
    simon_trl_idx = 1:150;
    if strcmp(SBJs{s},'PFC04')
        simon_trl_idx(72) = [];% from whrc neural variance
    elseif strcmp(SBJs{s},'PFC05')
        simon_trl_idx(76) = [];% from whrempty
    elseif strcmp(SBJs{s},'PFC01')
        simon_trl_idx(26) = [];% from whrtrl
    end
    all_bad_ix = unique([man_trl_rej_ix{s}'; sbj_data.bhv.empty_ix; sbj_data.bhv.bad_resp_ix; sbj_data.bhv.bad_rt_ix]);
    simon_bad_ix = [];
    for t = 1:length(all_bad_ix)
        simon_bad_ix = [simon_bad_ix; find(simon_trl_idx==all_bad_ix(t))];
    end
    if isempty(simon_bad_ix)
        fprintf('%s: All bad trials already tossed!\n',SBJs{s});
    else
        fprintf(2,'%s: Removing %d trials!\n',SBJs{s},length(simon_bad_ix));
        % Remove from behavior
        bhv_fields = fieldnames(bhv);
        for f = 1:length(bhv_fields)
            if length(bhv.(bhv_fields{f}))==length(simon_trl_idx) && ~contains(bhv_fields{f},'_ix')
                bhv.(bhv_fields{f})(simon_bad_ix) = [];
            end
        end
        % Remove from neural
        tfr.powspctrm(simon_bad_ix,:,:,:) = [];
    end
end

end