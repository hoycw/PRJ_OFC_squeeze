%% Plot group and subject level power distributions to inform modeling %%
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/');
addpath('/Users/colinhoy/Code/PRJ_OFC_squeeze/scripts/utils/');
close all
clear all

%%
% Baseline/ITI:
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'Sn8t0_bhvz_nrlz_out4';
% Stimulus decision phase:
an_id = 'TFRmth_S1t2_madS8t0_f2t40'; stat_id = 'S5t15_bhvz_nrl0_out3';%'S5t15_bhvz_nrl0_out2_bt1k';%
% an_id = 'TFRmth_S1t2_madA8t1_f2t40'; stat_id = 'S5t15_bhvz_nrlz_out4';
% Pre-decision:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'Dn5t0_bhvz_nrlz_out4';
% Post-decision/feedback:
% an_id = 'TFRmth_D1t1_madS8t0_f2t40'; stat_id = 'D0t1_bhvz_nrlz_out4';% stat_id = 'D0t5_bhvz_nrlz_out4';%

pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);
eval(['run ' prj_dir 'scripts/stat_vars/' stat_id '_vars.m']);
if ~isfield(st,'nboots'); error('only run this for permutation bootstrap stat_ids!'); end

fig_dir   = [prj_dir 'results/TFR/' an_id '/LMM/' stat_id '/power_distributions/'];
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

%% Load data
bhvs       = cell(size(SBJs));
mdls       = cell(size(SBJs));
for s = 1:length(SBJs)
    % Load behavior
    load([prj_dir 'data/' SBJs{s} '/' SBJs{s} '_stim_preproc.mat'],'sbj_data');
    bhvs{s} = sbj_data.bhv;
    mdls{s} = sbj_data.mdl;
end

%% Load group model table
table_all_fname = [prj_dir 'data/GRP/GRP_' an_id '_' stat_id '_full_table_all.csv'];
fprintf('\tLoading %s...\n',table_all_fname);
table_all = readtable(table_all_fname);

%% Test for normality of distribution
for f = 1:length(pow_vars)
%     [lil(f),lil_pval(f),lil_kstat(f),lil_critval(f)] = lillietest(table_all.(pow_vars{f}));
    fprintf('%s Lilliefors test: stat = %.3f, p = %.4f\n',pow_vars{f},lil_kstat(f),lil_pval(f));
end

%% Toss outliers
out_idx_all = struct;
out_ix_all = [];
good_tbl_all = struct;
for f = 1:length(pow_vars)
    % Identify outliers
    if isfield(st,'outlier_grp') && strcmp(st.outlier_grp,'sbj')
        for s = 1:length(SBJs)
            sbj_pow = table_all.(pow_vars{f})(sbj_idx);
            out_idx_all.(pow_vars{f})(table_all.sbj_n==s) = abs(sbj_pow) > ...
                abs(mean(sbj_pow))+(st.outlier_thresh*std(sbj_pow));
        end
    else
        out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>mean(table_all.(pow_vars{f}))...
            +(st.outlier_thresh*std(table_all.(pow_vars{f})));
    end
    
    % Toss outlier trials for each ROI and frequency band
    good_tbl_all.(pow_vars{f}) = table_all(~out_idx_all.(pow_vars{f}),:);
    
    % Report results
    fprintf('\n ================== Trials tossed for %s ==================\n',pow_vars{f});
    if any(out_idx_all.(pow_vars{f}))
        out_ix_all = [out_ix_all; find(out_idx_all.(pow_vars{f}))];
        fprintf(2,'\t%d outliers for %s in table_all:\t',sum(out_idx_all.(pow_vars{f})),pow_vars{f});
        fprintf(2,'%.2f, ',table_all.(pow_vars{f})(out_idx_all.(pow_vars{f})));
        fprintf('\n');
    else
        fprintf('No bad trials for %s with threshold %d\n',pow_vars{f},st.outlier_thresh);
    end
    fprintf('\tgood vs. all trials for %s in table_all = %d / %d\n',...
        pow_vars{f},size(good_tbl_all.(pow_vars{f}),1),size(table_all,1));
end
all_outliers_all = unique(out_ix_all);
fprintf(2,'Total bad trials in table_all: %d\n',length(all_outliers_all));

%% Plot distributions and outlier identification
pow_vars = {'PFC_theta','PFC_betalo','BG_theta','BG_betalo'};
sd_threshs = 2:4; out_colors = {'r','g','k'};
n_bins = 50;
for f = 1:length(pow_vars)
    fig_name = ['GRP_' an_id '_' stat_id '_' pow_vars{f} '_distributions'];
    figure('Name',fig_name,'units','norm','outerposition',[0 0 0.8 0.8]);
    n_out  = nan([length(sd_threshs) length(SBJs)+1]);
    thresh = nan([length(sd_threshs) length(SBJs)+1]);
    bins = linspace(min(table_all.(pow_vars{f}))-2, max(table_all.(pow_vars{f}))+2, n_bins);
    for s = 1:length(SBJs)+1
        if s>length(SBJs)
            pow_data = table_all.(pow_vars{f}); lab = 'GRP';
        else
            pow_data = table_all.(pow_vars{f})(table_all.sbj_n==s); lab = SBJs{s};
        end
        subplot(2,3,s); hold on;
        histogram(pow_data,bins);
        out_lines = []; out_lab = {};
        for o = 1:length(sd_threshs)
            if mean(pow_data)<0; error('negative mean!'); end
            thresh(o,s) = mean(pow_data)+(sd_threshs(o)*std(pow_data));
            n_out(o,s) = sum(abs(pow_data)>thresh(o,s));
            out_lines(o) = line([thresh(o,s) thresh(o,s)],ylim,'Color',out_colors{o});
            out_lab{o} = [num2str(sd_threshs(o)) 'SD=' num2str(thresh(o,s),'%.1f') ' (' num2str(n_out(o,s)) ' out)'];
        end
        xlim([min(table_all.(pow_vars{f}))-2 max(table_all.(pow_vars{f}))+2]);
        xlabel([pow_vars{f} ' (z)'],'Interpreter','none');
        title(sprintf('%s (mn=%.1f,sd=%.1f,n=%d)',lab,mean(pow_data),std(pow_data),length(pow_data)));
        set(gca,'FontSize',16);
        legend(out_lines,out_lab,'location','northeast');
    end
    subplot(2,3,6); hold on;
    sbj_lines = [];
    for s = 1:length(SBJs)
        sbj_lines(s) = line(thresh(:,s),n_out(:,s),'Color',sbj_colors(s,:));
    end
    sbj_lines(length(SBJs)+1) = line(thresh(:,end),n_out(:,end),'Color','k');
    legend(sbj_lines,[SBJs {'GRP'}],'location','best');
    xlabel('Thresholds');
    ylabel('# rejected trials');
    title(sprintf('sum of SBJ-level out=%d,%d,%d',sum(n_out(:,1:length(SBJs)),2)));
    set(gca,'FontSize',16);
    
    if save_fig
        fig_name = get(gcf,'Name');
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end


