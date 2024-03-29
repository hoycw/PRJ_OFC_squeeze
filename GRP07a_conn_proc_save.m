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
conn_metric = 'coh';%'PLVft';%'coh';%'PLV';%'ampcorr';
% an_ids = {'TFRmth_S03t2_f2t30_fourier'};%'TFRmth_S1t2_madS8t0_f2t40'};%'TFRmth_S1t2_madA8t1_f2t40'};%
an_ids = {'TFRmth_D1t1_f2t40_fourier'};%'TFRmth_D1t1_madS8t0_f2t40'};

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
        try
            fprintf('Loading %s\n',tfr_fname);
            load(tfr_fname,'tfr');
        catch
            error(['TFR for ampcorr does not exist, run GRP01a for ' an_id ' in ' SBJs{s}]);
        end
        if ~strcmp(conn_metric,'ampcorr') && ~strcmp(cfg_tfr.output,'fourier')
            error(['TFR output should be fourier for metric ' conn_metric]);
        end
        
        % Compute connectivity
        if contains(conn_metric,'coh')
            cfg = [];
            cfg.channelcmb = tfr.label';
            cfg.method = 'coh';
            conn = ft_connectivityanalysis(cfg,tfr);
        else
            conn = fn_connectivity_TFR(tfr,conn_metric);
        end
        
        % Compute jackknife connectivity
        if strcmp(conn_metric,'cohjk')
            % Jack-knife coherence for each trial contribution
            %   Recompute coherence after excluding one trial, then subtract
            %   from coherence with all trials to get single-trial contribution
            tic;
            conn_jk = nan([size(tfr.fourierspctrm,1) size(conn.cohspctrm)]);
            fprintf('Computing jackknife coherence (%d trials): ',size(tfr.fourierspctrm,1));
            for trl_ix = 1:size(tfr.fourierspctrm,1)
                if mod(trl_ix,5)==1; fprintf('\n\ntrial %d (%.2f min)...\n',trl_ix,toc/60); end
                jk_idx = setdiff(1:size(tfr.fourierspctrm,1),trl_ix);
                cfg.trials = jk_idx;
                tmp_coh = ft_connectivityanalysis(cfg,tfr);
                conn_jk(trl_ix,:,:,:,:) = conn.cohspctrm - tmp_coh.cohspctrm;
            end
            fprintf('\n');
            conn.cohspctrm = conn_jk;
            conn.dimord = ['rpt_' conn.dimord];
        end
        
        % SAVE data
        conn_fname = [sbj_dir SBJs{s} '_' an_ids{an_ix} '_' conn_metric '.mat'];
        fprintf('Saving %s\n',conn_fname);
        save(conn_fname,'-v7.3','conn');
        
        clear an an_vars_cmd conn tfr tfr_raw sbj_data conn conn_jk
    end
end
