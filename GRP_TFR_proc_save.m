%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
% clear all
% close all
% clc

prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
addpath([prj_dir 'scripts/']);
addpath([prj_dir 'scripts/utils/']);
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
SBJs = {'PFC03','PFC04','PFC05','PFC01'}; % 'PMC10'
% sbj_pfc_roi  = {'FPC', 'OFC', 'OFC', 'FPC'};

an_ids     = {'TFRw_S25t201_zbt25t05_fl2t40_varWid','simon_zbt','TFRw_S25t201_zbt25t05_fl2t40'};%,'TFRm_S25t201_zbtS_sm0_l0_wnVar','simon'};
% an_ids = {'TFRw_D101t201_zbt25t05_fl2t40_varWid','simon_D101t201_zbt25t05'};
% bsln_type  = 'zboot';
% bsln_boots = 500;

%% Time Frequency analysis
for s = 1:4
    %% Load data
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    tfr_fname = [sbj_dir SBJs{s} '_stim_preproc.mat']; % This is [-3 10] to cover all possible times
    fprintf('Loading %s\n',tfr_fname);
    load(tfr_fname,'sbj_data');
    
    %% Time-frequency representation
    for an_ix = 1:length(an_ids)
        an_vars_cmd = ['run ' prj_dir '/scripts/an_vars/' an_ids{an_ix} '_vars.m'];
        eval(an_vars_cmd);
        
%         % Trim to analysis time
%         if strcmp(an.event_type,'S')
%             [trial_lim_s_pad] = fn_get_filter_padding(cfg_tfr,an.trial_lim_s,an.bsln_lim);
%             cfgs = []; cfgs.latency = trial_lim_s_pad;
%             trim = ft_selectdata(cfgs,sbj_data.ts);
%         else
%             % That function doesn't work with variable length trials
%             trim = sbj_data.ts;
%         end
        
        % Extract time-frequency representation
        tfr_raw = ft_freqanalysis(cfg_tfr, sbj_data.ts);%trim);
        
        % Baseline correction
        switch an.bsln_type
            case {'zscore', 'zboot', 'demean', 'my_relchange'}
                tfr_bsln = fn_bsln_ft_tfr(tfr_raw,an.bsln_lim,an.bsln_type,an.bsln_boots);
            case {'relchange','db'}
                cfgbsln = [];
                cfgbsln.baseline     = an.bsln_lim;
                cfgbsln.baselinetype = an.bsln_type;
                cfgbsln.parameter    = 'powspctrm';
                tfr_bsln = ft_freqbaseline(cfgbsln,tfr_raw);
            case 'none'
                if an.complex
                    fprintf('\tNo baseline correction for ITPC data...\n');
                else
                    error('Why no baseline correction if not ITPC?');
                end
            otherwise
                error(['No baseline implemented for bsln_type: ' an.bsln_type]);
        end
        
        % Re-align to decision and trim
        if ~strcmp(an.event_type,an.bsln_evnt)
            if strcmp(an.event_type,'D') && strcmp(an.bsln_evnt,'S')
                [tfr_bsln] = fn_realign_tfr_s2r(tfr_bsln,sbj_data.bhv.rt,an.trial_lim_s);
            else
                error('unknown combination of analysis and baseline events');
            end
        end
        
        % Average trials and trim back down to analysis window
        cfgs = [];
        % cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
        cfgs.channel = an.ROI;
        cfgs.latency = an.trial_lim_s;
        tfr = ft_selectdata(cfgs, tfr_bsln);
        
        % SAVE data
        tfr_fname = [sbj_dir SBJs{s} '_' an_ids{an_ix} '.mat'];
        fprintf('Saving %s\n',tfr_fname);
        save(tfr_fname,'-v7.3','tfr');
    end
    
end