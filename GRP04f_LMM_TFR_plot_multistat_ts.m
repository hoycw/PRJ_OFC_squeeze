close all
clear all

%%
an_ids = {'TFRmth_S1t2_madA8t1_f2t40','TFRmth_S1t2_madS8t0_f2t40','TFRmth_D1t1_madS8t0_f2t40','TFRmth_D1t1_madS8t0_f2t40'};
stat_ids = {'Sn8t0_bhvz_nrlz_out4','S5t15_bhvz_nrlz_out4','Dn5t0_bhvz_nrlz_out4','D0t5_bhvz_nrlz_out4'};
if length(an_ids)~=length(stat_ids); error('an_ids and stat_ids must match!'); end

n_quantiles = 5;
save_fig = 1;
fig_ftype = 'png';

% Load SBJs, sbj_pfc_roi, sbj_bg_roi, and sbj_colors:
prj_dir = '/Users/colinhoy/Code/PRJ_OFC_squeeze/';
eval(['run ' prj_dir 'scripts/SBJ_vars.m']);


%% Compute LMMs
coef_ts = nan(size(stat_ids));
pval_ts = nan(size(stat_ids));
coef_fts = nan(size(stat_ids));
pval_fts = nan(size(stat_ids));
aic_fts = nan(size(stat_ids));
aic_ts = nan(size(stat_ids));
win_time = nan(size(stat_ids));
eta2f = nan(size(stat_ids));
rewp_eta2f = nan(size(stat_ids));
rewc_eta2f = nan(size(stat_ids));
% rewp_eta2 = nan(size(stat_ids));
for st_ix = 1:length(stat_ids)
    %% Load group model table
    eval(['run ' prj_dir 'scripts/stat_vars/' stat_ids{st_ix} '_vars.m']);
    table_all_fname = [prj_dir 'data/GRP/GRP_' an_ids{st_ix} '_' stat_ids{st_ix} '_full_table_all.csv'];
    fprintf('\tLoading %s...\n',table_all_fname);
    table_all = readtable(table_all_fname);
    win_time(st_ix) = mean(st.stat_lim);
    
    %% Toss outliers
    pow_vars = {'PFC_theta'};
    out_idx_all = struct;
    out_ix_all = [];
    good_tbl_all = struct;
    for f = 1:length(pow_vars)
        % Identify outliers
        out_idx_all.(pow_vars{f}) = abs(table_all.(pow_vars{f}))>st.outlier_thresh;
        
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
    
    %% Create previous trial and GRS tables
    % Toss NaNs from previous table
    good_tbl_prv = good_tbl_all;
    tbl_fields = good_tbl_all.(pow_vars{1}).Properties.VariableNames;
    for p = 1:length(pow_vars)
        prv_nan_idx = isnan(good_tbl_prv.(pow_vars{p}).SV_prv);
        good_tbl_prv.(pow_vars{p})(prv_nan_idx,:) = [];
        for f = 1:length(tbl_fields)
            if any(isnan(good_tbl_prv.(pow_vars{p}).(tbl_fields{f}))) && ~strcmp(tbl_fields{f},'grs')
                error(['NaN is table_prv.' tbl_fields{f}]);
            end
        end
    end
    
    %% Full model no decision ease/difficulty
    lme_full = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + reward_prv + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
    lme_full_norewp = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_cur + effortS_cur + effortS_prv + (1|sbj_n) + (1|trl_n_cur)');
    
    pfc_theta_rewp = compare(lme_full_norewp,lme_full,'CheckNesting',true)
    
    coef_fts(st_ix) = lme_full.Coefficients.Estimate(4);
    pval_fts(st_ix) = pfc_theta_rewp.pValue(2);
    aic_fts(st_ix) = lme_full.ModelCriterion.AIC;
    
    %% Compute effect size
    beta = fixedEffects(lme_full);
    y_est= [good_tbl_prv.PFC_theta.reward_cur good_tbl_prv.PFC_theta.effortS_cur ...
            good_tbl_prv.PFC_theta.reward_prv good_tbl_prv.PFC_theta.effortS_prv]*beta(2:5);
    y_rewp = good_tbl_prv.PFC_theta.reward_prv*beta(4);
    y_rewc = good_tbl_prv.PFC_theta.reward_cur*beta(2);
    residual = good_tbl_prv.PFC_theta.PFC_theta - y_est;
    eta2f(st_ix) = var(y_est)/var(good_tbl_prv.PFC_theta.PFC_theta); % total effect size
    rewp_eta2f(st_ix) = var(y_rewp)/(var(y_rewp)+var(residual)); % partial effect size for effect of x2 on y
    rewc_eta2f(st_ix) = var(y_rewc)/(var(y_rewc)+var(residual)); % partial effect size for effect of x2 on y
    
    %% PFC theta and previous reward:
    lme0 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ 1 + (1|sbj_n)');%,'StartMethod','random');
    lme1 = fitlme(good_tbl_prv.PFC_theta,'PFC_theta~ reward_prv + (1|sbj_n)');%,'StartMethod','random');
    pfc_theta_rew_prv = compare(lme0,lme1,'CheckNesting',true)%,'NSim',1000)
    
    coef_ts(st_ix) = lme1.Coefficients.Estimate(2);
    pval_ts(st_ix) = pfc_theta_rew_prv.pValue(2);
    aic_ts(st_ix) = lme1.ModelCriterion.AIC;
end

%% Plot the results and model fits effect sizes
% fig_dir   = [prj_dir 'results/TFR/' an_ids{ '/LMM/' stat_id '/PFC_theta/'];
% if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
figure;
subplot(4,1,1); hold on;
full_coef_line = plot(1:4,coef_fts,'r');
single_coef_line = plot(1:4,coef_ts,'b');
xlim([0 5]);
xticks([1:4]);
xticklabels({'Baseline','Stim','Pre-Dec','Post-Dec'});
ylabel('coefficient');
legend([full_coef_line,single_coef_line],{'full','single'},'Location','best');
set(gca,'FontSize',16);

subplot(4,1,2); hold on;
full_aic_line = plot(1:4,aic_fts,'r');
single_aic_line = plot(1:4,aic_ts,'b');
xlim([0 5]);
xticks([1:4]);
xticklabels({'Baseline','Stim','Pre-Dec','Post-Dec'});
% title('S-locked');
ylabel('AIC');
set(gca,'FontSize',16);
legend([full_aic_line,single_aic_line],{'full','single'},'Location','best');

subplot(4,1,3); hold on;
full_eta2_line = plot(1:4,eta2f,'r');
single_eta2_line = plot(1:4,eta2f,'b');
xlim([0 5]);
xticks([1:4]);
xticklabels({'Baseline','Stim','Pre-Dec','Post-Dec'});
% title('S-locked');
ylabel('Eta2');
set(gca,'FontSize',16);
legend([full_eta2_line,single_eta2_line],{'full','single'},'Location','best');

subplot(4,1,4); hold on;
full_rpeta2_line = plot(1:4,rewp_eta2f,'r');
single_rpeta2_line = plot(1:4,rewp_eta2f,'b');
full_rceta2_line = plot(1:4,rewc_eta2f,'m');
single_rceta2_line = plot(1:4,rewc_eta2f,'c');
xlim([0 5]);
xticks([1:4]);
xticklabels({'Baseline','Stim','Pre-Dec','Post-Dec'});
% title('S-locked');
ylabel('Partial Eta2');
set(gca,'FontSize',16);
legend([full_rpeta2_line,single_rpeta2_line,full_rceta2_line,single_rceta2_line],...
    {'full RewP','RewP','Full RewC','RewC'},'Location','best');
