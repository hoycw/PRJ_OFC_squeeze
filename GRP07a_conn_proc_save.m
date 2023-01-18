%% Colin preprocessing script
% Based on /Analysis/Group/Preprocess/Squeeze_preprocess_all_data_stimulus_locked.m
clear all
close all
% clc

addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
addpath('/Users/colinhoy/Code/Apps/fieldtrip/');
ft_defaults

%% Parameters
conn_metric = 'ampcorr';
an_ids = {'TFRmth_S1t2_madA8t1_f2t40'};%'TFRmth_S1t2_madS8t0_f2t40'};
% an_ids = {'TFRmth_D1t1_madS8t0_f2t40'};

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);

%% Time Frequency analysis
for s = 1:4
    sbj_dir = [prj_dir 'data/' SBJs{s} '/'];
    for an_ix = 1:length(an_ids)
        % Load TFR
        an_vars_cmd = ['run ' prj_dir 'scripts/an_vars/' an_ids{an_ix} '_vars.m'];
        eval(an_vars_cmd);
        
        tfr_fname = [sbj_dir SBJs{s} '_' an_ids{an_ix} '.mat'];
        if strcmp(conn_metric,'ampcorr')
            try 
                fprintf('Loading %s\n',tfr_fname);
                load(tfr_fname,'tfr');
            catch
                error(['TFR for ampcorr does not exist, run GRP01a for ' an_id]);
            end
        else
            if ~strcmp(cfg_tfr.output,'fourier')
                error(['TFR output should be fourier for metric ' conn_metric]);
            end
            % Load data
            preproc_fname = [sbj_dir SBJs{s} '_stim_preproc.mat']; % This is [-3 10] to cover all possible times
            fprintf('Loading %s\n',preproc_fname);
            load(preproc_fname,'sbj_data');
            
            % Compute complex TFR
            tfr_raw = ft_freqanalysis(cfg_tfr, sbj_data.ts);
            
            % Re-align to decision and re-trim
            if strcmp(an.event_type,'S')
                cfg_trim = [];
                cfg_trim.latency = an.trial_lim_s;
                tfr = ft_selectdata(cfg_trim,tfr_raw);
            elseif strcmp(an.event_type,'D')
                [tfr] = fn_realign_tfr_s2r(tfr_raw,sbj_data.bhv.rt,an.trial_lim_s);
            else
                error('only prepared for S or D-locked');
            end
        end
        
        % Compute connectivity
        conn = fn_connectivity_TFR(tfr,conn_metric);
        
        % SAVE data
        conn_fname = [sbj_dir SBJs{s} '_' an_ids{an_ix} '_' conn_metric '.mat'];
        fprintf('Saving %s\n',conn_fname);
        save(conn_fname,'-v7.3','conn');
        
        clear an an_vars_cmd conn tfr tfr_raw sbj_data
    end
end
