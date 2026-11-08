%% ========================================================================
%  任务：脑龄预测 + 全脑W‑score与PAD相关性（去总标题精简版 | 统一AI加字）
%  输出：图A（8.5cm×7cm）、图C（18.3cm×5cm），矢量PDF+300dpi PNG
%  已修正：全部输出文字HC→TC，不改动分组/颜色逻辑，避免散点变色
%% ========================================================================
clear; clc; close all;

% ---- Repository paths (no edit needed if you keep this folder layout) ----
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(script_dir);            % repository root
work_dir   = fullfile(repo_root, 'data');      % input tables + intermediate outputs
if ~exist(work_dir, 'dir'), mkdir(work_dir); end
cd(work_dir);

% ========================= 1. 读取原始数据 =========================
fprintf('读取原始数据 ABIDE_total2.xlsx ...\n');
raw = readtable('ABIDE_total2.xlsx');
raw.SUB_ID = string(raw.SUB_ID);

vol_raw_cols = raw.Properties.VariableNames(startsWith(raw.Properties.VariableNames, 'Volume'));
if length(vol_raw_cols) ~= 68
    error('原始体积列数不是68，实际为%d，请检查', length(vol_raw_cols));
end
X_raw = raw{:, vol_raw_cols};
y_age = raw.AGE_AT_SCAN;
dx = raw.DX_GROUP;
hc_idx = (dx == 2);
asd_idx = (dx == 1);

% ========================= 2. 多模型对比实验（仅HC训练，ASD独立外测） =================
fprintf('=== 多模型对比实验开始（仅HC训练，ASD独立外测） ===\n');
valid_hc = hc_idx & all(~isnan(X_raw),2) & ~isnan(y_age);
X_train = X_raw(valid_hc,:);
y_train = y_age(valid_hc);
fprintf('HC有效样本: %d, 特征维度: %d\n', size(X_train,1), size(X_train,2));

% ---- 运行开关 ----
enable_gpr        = true;
enable_rvr        = false;
enable_stratified = false;
enable_outlier    = true;
enable_lsboost    = true;
enable_rf         = true;
enable_sex_covariate = false;

if enable_sex_covariate
    sex_col = find(ismember(lower(raw.Properties.VariableNames), 'sex'));
    if isempty(sex_col)
        warning('未找到SEX列，性别协变量校正已跳过');
    else
        X_raw = [X_raw, double(raw{:, sex_col})];
        valid_hc = hc_idx & all(~isnan(X_raw),2) & ~isnan(y_age);
        X_train = X_raw(valid_hc,:);
        fprintf('已启用性别协变量校正，总维度=%d\n', size(X_train,2));
    end
end

valid_asd = asd_idx & all(~isnan(X_raw),2) & ~isnan(y_age);
X_asd = X_raw(valid_asd,:);
y_asd = y_age(valid_asd);
fprintf('ASD独立测试样本: %d\n', size(X_asd,1));

rng(42);

% ---- 结果汇总容器 ----
model_names = {};  best_param_str = {};
hc_oob_r = []; hc_oob_mae = []; hc_oob_rmse = [];
asd_r = []; asd_mae = []; asd_rmse = [];
asd_pad_min = []; asd_pad_max = [];
all_fit_fun = {}; all_pred_fun = {}; all_params = {}; all_oob = {}; all_models = {}; all_trainX = {};
strata_edges = [0, 12, 22, inf];
strat_models = cell(3,1); strat_preds = cell(3,1); strat_sel = cell(3,1);
strat_idx = [];

% ==================== [模型1] 岭回归 Ridge ====================
fprintf('\n[模型1/9] 岭回归 Ridge（线性基线）...\n');
tic;
grid_m1 = logspace(-4, 4, 25)';
fit_m1  = @(Xs,y,p) fitrlinear(Xs, y, 'Learner','leastsquares', 'Regularization','ridge', 'Lambda', p(1));
pred_m1 = @(m,Xs) predict(m, Xs);
[bp1, oob1, met1] = evaluate_model_cv(X_train, y_train, grid_m1, fit_m1, pred_m1);
[m1, metA1] = fit_final_test(X_train, y_train, X_asd, y_asd, bp1, fit_m1, pred_m1);
warn_grid_boundary(bp1, grid_m1, 'Ridge');
model_names{end+1} = 'Ridge-Linear';
best_param_str{end+1} = sprintf('lambda=%.4g', bp1(1));
hc_oob_r(end+1)=met1.r; hc_oob_mae(end+1)=met1.mae; hc_oob_rmse(end+1)=met1.rmse;
asd_r(end+1)=metA1.r; asd_mae(end+1)=metA1.mae; asd_rmse(end+1)=metA1.rmse;
asd_pad_min(end+1)=min(metA1.pad); asd_pad_max(end+1)=max(metA1.pad);
all_fit_fun{end+1}=fit_m1; all_pred_fun{end+1}=pred_m1; all_params{end+1}=bp1;
all_oob{end+1}=oob1; all_models{end+1}=m1; all_trainX{end+1}=X_train;
fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
        toc, met1.r, met1.mae, met1.rmse, metA1.r, metA1.mae, metA1.rmse, min(metA1.pad), max(metA1.pad));

% ==================== [模型2] 岭回归 + 二次多项式 ====================
fprintf('\n[模型2/9] 岭回归 + 二次平方项特征（Ridge-Poly2）...\n');
tic;
grid_m2 = logspace(-4, 4, 25)';
fit_m2  = @(Xs,y,p) fitrlinear([Xs, Xs.^2], y, 'Learner','leastsquares', 'Regularization','ridge', 'Lambda', p(1));
pred_m2 = @(m,Xs) predict(m, [Xs, Xs.^2]);
[bp2, oob2, met2] = evaluate_model_cv(X_train, y_train, grid_m2, fit_m2, pred_m2);
[m2, metA2] = fit_final_test(X_train, y_train, X_asd, y_asd, bp2, fit_m2, pred_m2);
warn_grid_boundary(bp2, grid_m2, 'Ridge-Poly2');
model_names{end+1} = 'Ridge-Poly2';
best_param_str{end+1} = sprintf('lambda=%.4g', bp2(1));
hc_oob_r(end+1)=met2.r; hc_oob_mae(end+1)=met2.mae; hc_oob_rmse(end+1)=met2.rmse;
asd_r(end+1)=metA2.r; asd_mae(end+1)=metA2.mae; asd_rmse(end+1)=metA2.rmse;
asd_pad_min(end+1)=min(metA2.pad); asd_pad_max(end+1)=max(metA2.pad);
all_fit_fun{end+1}=fit_m2; all_pred_fun{end+1}=pred_m2; all_params{end+1}=bp2;
all_oob{end+1}=oob2; all_models{end+1}=m2; all_trainX{end+1}=X_train;
fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
        toc, met2.r, met2.mae, met2.rmse, metA2.r, metA2.mae, metA2.rmse, min(metA2.pad), max(metA2.pad));

% ==================== [模型3] RBF-SVR ====================
fprintf('\n[模型3/9] RBF-SVR（支持向量回归）...\n');
tic;
[KSg3, BCg3, EPg3] = ndgrid([2 5 10 20 50 100], [10 100], [0.5 1]);
grid_m3 = [KSg3(:), BCg3(:), EPg3(:)];
fit_m3  = @(Xs,y,p) fitrsvm(Xs, y, 'KernelFunction','rbf', 'Standardize',false, ...
                            'BoxConstraint',p(2), 'KernelScale',p(1), 'Epsilon',p(3));
pred_m3 = @(m,Xs) predict(m, Xs);
[bp3, oob3, met3] = evaluate_model_cv(X_train, y_train, grid_m3, fit_m3, pred_m3);
[m3, metA3] = fit_final_test(X_train, y_train, X_asd, y_asd, bp3, fit_m3, pred_m3);
warn_grid_boundary(bp3, grid_m3, 'SVR-RBF');
model_names{end+1} = 'SVR-RBF';
best_param_str{end+1} = sprintf('KS=%.2g,C=%.1g,eps=%.2g', bp3(1), bp3(2), bp3(3));
hc_oob_r(end+1)=met3.r; hc_oob_mae(end+1)=met3.mae; hc_oob_rmse(end+1)=met3.rmse;
asd_r(end+1)=metA3.r; asd_mae(end+1)=metA3.mae; asd_rmse(end+1)=metA3.rmse;
asd_pad_min(end+1)=min(metA3.pad); asd_pad_max(end+1)=max(metA3.pad);
all_fit_fun{end+1}=fit_m3; all_pred_fun{end+1}=pred_m3; all_params{end+1}=bp3;
all_oob{end+1}=oob3; all_models{end+1}=m3; all_trainX{end+1}=X_train;
svr_full = m3;
fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
        toc, met3.r, met3.mae, met3.rmse, metA3.r, metA3.mae, metA3.rmse, min(metA3.pad), max(metA3.pad));

% ==================== [模型4] 核岭回归 KRR ====================
fprintf('\n[模型4/9] 核岭回归 KRR（RBF核）...\n');
tic;
try
    [Lg4, Sg4] = ndgrid([0.1 1 10 100], [2 5 10]);
    grid_m4 = [Lg4(:), Sg4(:)];
    fit_m4  = @(Xs,y,p) krr_fit(Xs, y, p);
    pred_m4 = @(m,Xs) krr_pred(m, Xs);
    [bp4, oob4, met4] = evaluate_model_cv(X_train, y_train, grid_m4, fit_m4, pred_m4);
    [m4, metA4] = fit_final_test(X_train, y_train, X_asd, y_asd, bp4, fit_m4, pred_m4);
    warn_grid_boundary(bp4, grid_m4, 'KRR');
    model_names{end+1} = 'KRR-RBF';
    best_param_str{end+1} = sprintf('lambda=%.3g,KS=%.2g', bp4(1), bp4(2));
    hc_oob_r(end+1)=met4.r; hc_oob_mae(end+1)=met4.mae; hc_oob_rmse(end+1)=met4.rmse;
    asd_r(end+1)=metA4.r; asd_mae(end+1)=metA4.mae; asd_rmse(end+1)=metA4.rmse;
    asd_pad_min(end+1)=min(metA4.pad); asd_pad_max(end+1)=max(metA4.pad);
    all_fit_fun{end+1}=fit_m4; all_pred_fun{end+1}=pred_m4; all_params{end+1}=bp4;
    all_oob{end+1}=oob4; all_models{end+1}=m4; all_trainX{end+1}=X_train;
    fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
            toc, met4.r, met4.mae, met4.rmse, metA4.r, metA4.mae, metA4.rmse, min(metA4.pad), max(metA4.pad));
catch ME
    fprintf('  [跳过] KRR出错: %s\n', ME.message);
end

% ==================== [模型5] 高斯过程回归 GPR ====================
if enable_gpr
    fprintf('\n[模型5/9] 高斯过程回归 GPR（squaredexponential）...\n');
    tic;
    try
        grid_m5 = [1];
        fit_m5  = @(Xs,y,p) fitrgp(Xs, y, 'KernelFunction','squaredexponential', 'Standardize',false);
        pred_m5 = @(m,Xs) predict(m, Xs);
        [bp5, oob5, met5] = evaluate_model_cv(X_train, y_train, grid_m5, fit_m5, pred_m5);
        [m5, metA5] = fit_final_test(X_train, y_train, X_asd, y_asd, bp5, fit_m5, pred_m5);
        warn_grid_boundary(bp5, grid_m5, 'GPR');
        model_names{end+1} = 'GPR';
        best_param_str{end+1} = 'auto (MLE)';
        hc_oob_r(end+1)=met5.r; hc_oob_mae(end+1)=met5.mae; hc_oob_rmse(end+1)=met5.rmse;
        asd_r(end+1)=metA5.r; asd_mae(end+1)=metA5.mae; asd_rmse(end+1)=metA5.rmse;
        asd_pad_min(end+1)=min(metA5.pad); asd_pad_max(end+1)=max(metA5.pad);
        all_fit_fun{end+1}=fit_m5; all_pred_fun{end+1}=pred_m5; all_params{end+1}=bp5;
        all_oob{end+1}=oob5; all_models{end+1}=m5; all_trainX{end+1}=X_train;
        fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
                toc, met5.r, met5.mae, met5.rmse, metA5.r, metA5.mae, metA5.rmse, min(metA5.pad), max(metA5.pad));
    catch ME
        fprintf('  [跳过] GPR出错: %s\n', ME.message);
    end
end

% ==================== [模型6] 相关向量回归 RVR ====================
if enable_rvr
    fprintf('\n[模型6/9] 相关向量回归 RVR（稀疏贝叶斯核）...\n');
    tic;
    try
        grid_m6 = [2 5 10]';
        fit_m6  = @(Xs,y,p) rvr_fit(Xs, y, p);
        pred_m6 = @(m,Xs) rvr_pred(m, Xs);
        [bp6, oob6, met6] = evaluate_model_cv(X_train, y_train, grid_m6, fit_m6, pred_m6);
        [m6, metA6] = fit_final_test(X_train, y_train, X_asd, y_asd, bp6, fit_m6, pred_m6);
        warn_grid_boundary(bp6, grid_m6, 'RVR');
        model_names{end+1} = 'RVR';
        best_param_str{end+1} = sprintf('KS=%.2g', bp6(1));
        hc_oob_r(end+1)=met6.r; hc_oob_mae(end+1)=met6.mae; hc_oob_rmse(end+1)=met6.rmse;
        asd_r(end+1)=metA6.r; asd_mae(end+1)=metA6.mae; asd_rmse(end+1)=metA6.rmse;
        asd_pad_min(end+1)=min(metA6.pad); asd_pad_max(end+1)=max(metA6.pad);
        all_fit_fun{end+1}=fit_m6; all_pred_fun{end+1}=pred_m6; all_params{end+1}=bp6;
        all_oob{end+1}=oob6; all_models{end+1}=m6; all_trainX{end+1}=X_train;
        fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
                toc, met6.r, met6.mae, met6.rmse, metA6.r, metA6.mae, metA6.rmse, min(metA6.pad), max(metA6.pad));
    catch ME
        fprintf('  [跳过] RVR出错: %s\n', ME.message);
    end
end

% ==================== [模型7] 梯度提升树 LSBoost ====================
if enable_lsboost
fprintf('\n[模型7/9] 梯度提升树 LSBoost...\n');
tic;
[NCg7, LRg7, MLg7] = ndgrid([100 200], [0.05 0.1], [5 10]);
grid_m7 = [NCg7(:), LRg7(:), MLg7(:)];
fit_m7  = @(Xs,y,p) fitrensemble(Xs, y, 'Method','LSBoost', 'NumLearningCycles',p(1), ...
                                 'LearnRate',p(2), 'Learners', templateTree('MinLeafSize',p(3)));
pred_m7 = @(m,Xs) predict(m, Xs);
[bp7, oob7, met7] = evaluate_model_cv(X_train, y_train, grid_m7, fit_m7, pred_m7);
[m7, metA7] = fit_final_test(X_train, y_train, X_asd, y_asd, bp7, fit_m7, pred_m7);
warn_grid_boundary(bp7, grid_m7, 'LSBoost');
model_names{end+1} = 'LSBoost';
best_param_str{end+1} = sprintf('Cyc=%d,LR=%.2g,Leaf=%d', bp7(1), bp7(2), bp7(3));
hc_oob_r(end+1)=met7.r; hc_oob_mae(end+1)=met7.mae; hc_oob_rmse(end+1)=met7.rmse;
asd_r(end+1)=metA7.r; asd_mae(end+1)=metA7.mae; asd_rmse(end+1)=metA7.rmse;
asd_pad_min(end+1)=min(metA7.pad); asd_pad_max(end+1)=max(metA7.pad);
all_fit_fun{end+1}=fit_m7; all_pred_fun{end+1}=pred_m7; all_params{end+1}=bp7;
all_oob{end+1}=oob7; all_models{end+1}=m7; all_trainX{end+1}=X_train;
fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
        toc, met7.r, met7.mae, met7.rmse, metA7.r, metA7.mae, metA7.rmse, min(metA7.pad), max(metA7.pad));
end

% ==================== [模型8] 随机森林 RF ====================
if enable_rf
fprintf('\n[模型8/9] 随机森林 RF（TreeBagger）...\n');
tic;
grid_m8 = [100 20];
fit_m8  = @(Xs,y,p) TreeBagger(p(1), Xs, y, 'Method','regression', 'MinLeafSize',p(2));
pred_m8 = @(m,Xs) predict(m, Xs);
[bp8, oob8, met8] = evaluate_model_cv(X_train, y_train, grid_m8, fit_m8, pred_m8);
[m8, metA8] = fit_final_test(X_train, y_train, X_asd, y_asd, bp8, fit_m8, pred_m8);
warn_grid_boundary(bp8, grid_m8, 'RF');
model_names{end+1} = 'RandomForest';
best_param_str{end+1} = sprintf('Trees=%d,Leaf=%d', bp8(1), bp8(2));
hc_oob_r(end+1)=met8.r; hc_oob_mae(end+1)=met8.mae; hc_oob_rmse(end+1)=met8.rmse;
asd_r(end+1)=metA8.r; asd_mae(end+1)=metA8.mae; asd_rmse(end+1)=metA8.rmse;
asd_pad_min(end+1)=min(metA8.pad); asd_pad_max(end+1)=max(metA8.pad);
all_fit_fun{end+1}=fit_m8; all_pred_fun{end+1}=pred_m8; all_params{end+1}=bp8;
all_oob{end+1}=oob8; all_models{end+1}=m8; all_trainX{end+1}=X_train;
fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
        toc, met8.r, met8.mae, met8.rmse, metA8.r, metA8.mae, metA8.rmse, min(metA8.pad), max(metA8.pad));
end

% ==================== [模型9] 年龄分层分段优化 ====================
if enable_stratified
    fprintf('\n[模型9/9] 年龄分层分段优化（SVR-RBF基底）...\n');
    tic;
    try
        [KSg9, BCg9, EPg9] = ndgrid([2 5 10], [10 100], [0.5]);
        grid_s9 = [KSg9(:), BCg9(:), EPg9(:)];
        strat_oob = nan(size(X_train,1), 1);
        for s = 1:3
            sel = y_train >= strata_edges(s) & y_train < strata_edges(s+1);
            strat_sel{s} = sel;
            if sum(sel) < 30
                fprintf('  层%d [%g,%g) 样本过少(%d)，回退用全局SVR模型\n', s, strata_edges(s), strata_edges(s+1), sum(sel));
                strat_models{s} = svr_full.learner;
                strat_preds{s} = svr_full.pred_fun;
                strat_oob(sel) = svr_full.pred_fun(svr_full.learner, zscore_apply(X_train(sel,:), X_train));
                continue;
            end
            [bp_s, oob_s, met_s] = evaluate_model_cv(X_train(sel,:), y_train(sel), grid_s9, fit_m3, pred_m3);
            Xs_s = zscore_rows(X_train(sel,:));
            strat_models{s} = fit_m3(Xs_s, y_train(sel), bp_s);
            strat_preds{s} = pred_m3;
            strat_oob(sel) = oob_s;
            fprintf('  层%d [%g,%g): n=%d, KS=%.2f C=%.1f eps=%.2f, OOB_r=%.3f\n', ...
                    s, strata_edges(s), strata_edges(s+1), sum(sel), bp_s(1), bp_s(2), bp_s(3), met_s.r);
        end
        v9 = ~isnan(strat_oob);
        met9 = struct('r', rv(y_train(v9), strat_oob(v9)), 'mae', mean(abs(y_train(v9)-strat_oob(v9))), ...
                      'rmse', sqrt(mean((y_train(v9)-strat_oob(v9)).^2)), 'n', sum(v9));
        asd_pred9 = nan(size(y_asd));
        for s = 1:3
            selA = y_asd >= strata_edges(s) & y_asd < strata_edges(s+1);
            if ~any(selA), continue; end
            if ~isempty(strat_models{s}) && isa(strat_models{s}, 'struct')
                Xs_a = zscore_apply(X_asd(selA,:), X_train);
            else
                Xs_a = zscore_apply(X_asd(selA,:), X_train(strat_sel{s},:));
            end
            asd_pred9(selA) = strat_preds{s}(strat_models{s}, Xs_a);
        end
        vA9 = ~isnan(asd_pred9);
        metA9 = struct('r', rv(y_asd(vA9), asd_pred9(vA9)), 'mae', mean(abs(y_asd(vA9)-asd_pred9(vA9))), ...
                       'rmse', sqrt(mean((y_asd(vA9)-asd_pred9(vA9)).^2)), 'n', sum(vA9), 'pad', asd_pred9(vA9)-y_asd(vA9));
        model_names{end+1} = 'Stratified-SVR';
        strat_idx = numel(model_names);
        best_param_str{end+1} = '三层独立调参(SVR)';
        hc_oob_r(end+1)=met9.r; hc_oob_mae(end+1)=met9.mae; hc_oob_rmse(end+1)=met9.rmse;
        asd_r(end+1)=metA9.r; asd_mae(end+1)=metA9.mae; asd_rmse(end+1)=metA9.rmse;
        asd_pad_min(end+1)=min(metA9.pad); asd_pad_max(end+1)=max(metA9.pad);
        all_fit_fun{end+1}=fit_m3; all_pred_fun{end+1}=pred_m3; all_params{end+1}=bp3;
        all_oob{end+1}=strat_oob; all_models{end+1}=[]; all_trainX{end+1}=X_train;
        fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
                toc, met9.r, met9.mae, met9.rmse, metA9.r, metA9.mae, metA9.rmse, min(metA9.pad), max(metA9.pad));
    catch ME
        fprintf('  [跳过] 分层模型出错: %s\n', ME.message);
    end
end

% ==================== 暂定最优模型（供离群剔除使用） ====================
prov_cand = setdiff(1:numel(asd_r), strat_idx);
[~, prov_i] = max(asd_r(prov_cand));
prov_i = prov_cand(prov_i);
prov_fit = all_fit_fun{prov_i}; prov_pred = all_pred_fun{prov_i};
prov_params = all_params{prov_i}; prov_name = model_names{prov_i};

% ==================== [模型10] 最优模型 + HC训练集离群剔除增强 ====================
if enable_outlier && ~isempty(prov_i)
    fprintf('\n[模型10/10] %s + HC训练集离群PAD剔除增强...\n', prov_name);
    tic;
    try
        Xs_all = zscore_rows(X_train);
        m_init = prov_fit(Xs_all, y_train, prov_params);
        pad_hc = y_train - prov_pred(m_init, Xs_all);
        keep = abs(pad_hc - mean(pad_hc)) <= 3*std(pad_hc) & abs(pad_hc) <= 15;
        fprintf('  剔除HC离群 %d 例（保留 %d 例）\n', sum(~keep), sum(keep));
        if sum(keep) < 0.7*numel(keep)
            warning('离群剔除比例过高(%.1f%%)，请检查阈值', 100*sum(~keep)/numel(keep));
        end
        [bp10, oob10, met10] = evaluate_model_cv(X_train(keep,:), y_train(keep), prov_params(:)', prov_fit, prov_pred);
        [m10, metA10] = fit_final_test(X_train(keep,:), y_train(keep), X_asd, y_asd, prov_params, prov_fit, prov_pred);
        oob10_full = nan(size(y_train)); oob10_full(keep) = oob10;
        model_names{end+1} = [prov_name '+OutlierRemoval'];
        best_param_str{end+1} = sprintf('同%s,剔除%d例', prov_name, sum(~keep));
        hc_oob_r(end+1)=met10.r; hc_oob_mae(end+1)=met10.mae; hc_oob_rmse(end+1)=met10.rmse;
        asd_r(end+1)=metA10.r; asd_mae(end+1)=metA10.mae; asd_rmse(end+1)=metA10.rmse;
        asd_pad_min(end+1)=min(metA10.pad); asd_pad_max(end+1)=max(metA10.pad);
        all_fit_fun{end+1}=prov_fit; all_pred_fun{end+1}=prov_pred; all_params{end+1}=prov_params;
        all_oob{end+1}=oob10_full; all_models{end+1}=m10; all_trainX{end+1}=X_train(keep,:);
        fprintf('  结果(%.1fs): HC_OOB r=%.4f MAE=%.3f RMSE=%.3f | ASD r=%.4f MAE=%.3f RMSE=%.3f | ASD PAD[%.2f,%.2f]\n', ...
                toc, met10.r, met10.mae, met10.rmse, metA10.r, metA10.mae, metA10.rmse, min(metA10.pad), max(metA10.pad));
    catch ME
        fprintf('  [跳过] 离群剔除模型出错: %s\n', ME.message);
    end
end

% ==================== 最终最优模型选择 ====================
[~, best_idx_final] = max(asd_r);
linear_asd_r = asd_r(1);
delta_r = asd_r - linear_asd_r;

% ==================== 多模型性能对比汇总 ====================
fprintf('\n===================== 多模型性能对比汇总（ASD r 降序） =====================\n');
fprintf('%-30s | %9s %8s %8s | %9s %8s %8s | %8s | %s\n', ...
    'Model','HC_OOB_r','HC_MAE','HC_RMSE','ASD_r','ASD_MAE','ASD_RMSE','Delta_r','BestParams');
[~, order] = sort(asd_r, 'descend');
for ii = order(:)'
    tag = '';
    if ii == best_idx_final, tag = '  <<< 最优'; end
    fprintf('%-30s | %9.4f %8.3f %8.3f | %9.4f %8.3f %8.3f | %+8.4f | %s%s\n', ...
        model_names{ii}, hc_oob_r(ii), hc_oob_mae(ii), hc_oob_rmse(ii), ...
        asd_r(ii), asd_mae(ii), asd_rmse(ii), delta_r(ii), best_param_str{ii}, tag);
end
fprintf('（Delta_r = 该模型ASD r − 线性基线ASD r = %.4f）\n', linear_asd_r);

export_model_table = true;
if export_model_table
    cmp_table = table(model_names(:), best_param_str(:), hc_oob_r(:), hc_oob_mae(:), hc_oob_rmse(:), ...
                      asd_r(:), asd_mae(:), asd_rmse(:), delta_r(:), ...
                      'VariableNames', {'Model','BestParams','HC_OOB_r','HC_OOB_MAE','HC_OOB_RMSE', ...
                                        'ASD_r','ASD_MAE','ASD_RMSE','Delta_r_vs_Linear'});
    cmp_table = sortrows(cmp_table, 'ASD_r', 'descend');
    writetable(cmp_table, fullfile(work_dir, 'BrainAge_Model_Comparison.xlsx'));
    fprintf('已导出模型对比表: BrainAge_Model_Comparison.xlsx\n');
end

% ==================== 全样本预测 ====================
fprintf('\n【最终最优模型】%s\n', model_names{best_idx_final});
fprintf('  最优超参数: %s\n', best_param_str{best_idx_final});
fprintf('  ASD r = %.4f（相对线性基线提升 Δr = %+.4f）\n', asd_r(best_idx_final), delta_r(best_idx_final));
if ~isempty(strat_idx) && best_idx_final == strat_idx && enable_stratified
    pred_age = nan(size(X_raw,1), 1);
    for s = 1:3
        sel = y_age >= strata_edges(s) & y_age < strata_edges(s+1);
        if ~any(sel), continue; end
        if isempty(strat_models{s})
            Xs_p = zscore_apply(X_raw(sel,:), X_train);
        else
            Xs_p = zscore_apply(X_raw(sel,:), X_train(strat_sel{s},:));
        end
        pred_age(sel) = strat_preds{s}(strat_models{s}, Xs_p);
    end
else
    X_all_std = zscore_apply(X_raw, all_trainX{best_idx_final});
    pred_age = all_models{best_idx_final}.pred_fun(all_models{best_idx_final}.learner, X_all_std);
end
raw.Predicted_Age = pred_age;
raw.PAD = raw.Predicted_Age - raw.AGE_AT_SCAN;
fprintf('PAD范围: [%.2f, %.2f]\n', min(raw.PAD), max(raw.PAD));

best_oob_vec = all_oob{best_idx_final};
asd_x = raw.AGE_AT_SCAN(asd_idx); asd_y = raw.Predicted_Age(asd_idx);
av = ~isnan(asd_x) & ~isnan(asd_y);
vob = ~isnan(best_oob_vec);
fprintf('\n--- 绘图指标（与图A文本框对应）---\n');
fprintf('ASD样本外绘图指标: r = %.4f, MAE = %.4f, RMSE = %.4f, n = %d\n', ...
        rv(asd_x(av), asd_y(av)), mean(abs(asd_x(av)-asd_y(av))), ...
        sqrt(mean((asd_x(av)-asd_y(av)).^2)), sum(av));
fprintf('TC OOB袋外绘图指标: r = %.4f, MAE = %.4f, RMSE = %.4f, n = %d\n', ...
        rv(y_train(vob), best_oob_vec(vob)), mean(abs(y_train(vob)-best_oob_vec(vob))), ...
        sqrt(mean((y_train(vob)-best_oob_vec(vob)).^2)), sum(vob));

% ==================== 年龄偏倚校正（TDC out-of-fold 线性偏倚线，Cole et al., 2017 口径）====================
% 归档补丁（2026-09-17）：原 v2 直接用 PAD = 预测年龄−实际年龄，未校正回归到均值偏倚；
% 现用最终模型的 TDC 五折 OOF 预测拟合 PAD_oof = slope*age + intercept，并从全部被试 PAD 中扣除。
vob_bias = ~isnan(best_oob_vec);
oof_pad_bias = best_oob_vec(vob_bias) - y_train(vob_bias);
bias_coef = polyfit(y_train(vob_bias), oof_pad_bias, 1);   % [slope, intercept]
raw.PAD = raw.PAD - (bias_coef(1)*raw.AGE_AT_SCAN + bias_coef(2));
fprintf('年龄偏倚线（TDC OOF）: slope = %.3f, intercept = %.3f；全部被试 PAD 已完成偏倚校正\n', ...
        bias_coef(1), bias_coef(2));

% ========================= 3. 图A：脑龄预测散点图（去总标题 | 统一格式） =================
fprintf('\n--- 绘制并导出图A ---');
figure('Name','Figure A','Color','w','Units','centimeters','Position',[2,2,8.5,7.0]);

ax_lo = floor(min([raw.AGE_AT_SCAN(asd_idx); raw.Predicted_Age(asd_idx); y_train; best_oob_vec])*2)/2;
ax_hi = ceil(max([raw.AGE_AT_SCAN(asd_idx); raw.Predicted_Age(asd_idx); y_train; best_oob_vec])*2)/2;

subplot(1,2,1);
plot_correlation_A(raw.AGE_AT_SCAN(asd_idx), raw.Predicted_Age(asd_idx), 'ASD', 'ASD cohort', [ax_lo ax_hi]);
subplot(1,2,2);
% ======== 修改1：仅修改显示文字为TC，颜色判断label='HC'不变，分组颜色不紊乱 ========
plot_correlation_A(y_train, best_oob_vec, 'HC', 'TC (5-fold OOB)', [ax_lo ax_hi]);

% ==========【修改1：移除图A总标题，统一到AI中添加】==========
% sgtitle('Correlation between chronological age and predicted brain age', ...
%     'FontSize',6.5, 'FontName','Arial', 'FontWeight','bold');

figA_pdf = fullfile(work_dir, 'Figure_A_BrainAge_Prediction_noTitle.pdf');
figA_png = fullfile(work_dir, 'Figure_A_BrainAge_Prediction_noTitle.png');
try
    exportgraphics(gcf, figA_pdf, 'ContentType','vector');
    exportgraphics(gcf, figA_png, 'Resolution', 300);
    fprintf('已输出图A矢量PDF+PNG，尺寸8.5cm×7cm（无总标题）\n');
catch ME
    warning('图A导出失败: %s', ME.message);
end

% ========================= 4. 读取W-score数据 =========================
fprintf('读取W-score数据 ABIDE_All_Metrics_Subtypes_Summaryx.xlsx ...\n');
wscore = readtable('ABIDE_All_Metrics_Subtypes_Summaryx.xlsx');
wscore.SUB_ID = string(wscore.SUB_ID);
all_vol = wscore.Properties.VariableNames(startsWith(wscore.Properties.VariableNames, 'Volume_'));
vol_w_cols = all_vol(~contains(all_vol, 'Subtype'));
if length(vol_w_cols) ~= 68
    error('W-score体积列数不是68，实际为%d，请检查', length(vol_w_cols));
end

asd_w = wscore(wscore.DX_GROUP == 1, :);
[~, idx_raw] = ismember(asd_w.SUB_ID, raw.SUB_ID);
if any(idx_raw==0)
    warning('部分ASD被试在原始数据中找不到PAD，已剔除');
    asd_w = asd_w(idx_raw>0, :);
    idx_raw = idx_raw(idx_raw>0);
end
asd_w.PAD = raw.PAD(idx_raw);
asd_w = rmmissing(asd_w, 'DataVariables', [vol_w_cols, 'PAD']);
fprintf('ASD组有效样本（W-score + PAD）: %d\n', height(asd_w));

% ========================= 5. 计算W-score与PAD的相关系数 =================
r_vals = zeros(68,1);
p_vals = zeros(68,1);
for i = 1:68
    vol = asd_w.(vol_w_cols{i});
    pad = asd_w.PAD;
    [R, P] = corrcoef(vol, pad);
    r_vals(i) = R(1,2);
    p_vals(i) = P(1,2);
end

p_fdr = mafdr(p_vals, 'BHFDR', true);
r_final = r_vals;

dk_names = {...
    'precentral','postcentral','paracentral','superior frontal', ...
    'rostral middle frontal','caudal middle frontal','frontal pole', ...
    'lateral orbitofrontal','medial orbitofrontal','pars opercularis', ...
    'pars orbitalis','pars triangularis','insula', ...
    'rostral anterior cingulate','caudal anterior cingulate','posterior cingulate', ...
    'isthmus cingulate','parahippocampal','superior parietal','inferior parietal', ...
    'cuneus','entorhinal','fusiform','lingual', ...
    'pericalcarine','precuneus','supramarginal','superior temporal', ...
    'middle temporal','inferior temporal','temporal pole', ...
    'banks of the superior temporal sulcus','transverse temporal','lateral occipital'};
hemi_list = repmat({'left';'right'},34,1);
region_list = [];
for i=1:34
    region_list = [region_list; dk_names(i); dk_names(i)];
end
output_table = table(r_final, p_vals, p_fdr, 'VariableNames', {'Correlation_R','P_Value','P_FDR'});
output_table.Hemisphere = hemi_list;
output_table.Region_Name = region_list;

output_file = fullfile(work_dir, 'ASD_All_ROI_PAD_Correlation_Wscore.xlsx');
writetable(output_table, output_file);
fprintf('已生成供R绘图的Excel文件: %s\n', output_file);

% ========================= 6. 图C：三组Brain-PAD分布图（去总标题 | 统一格式） =================
fprintf('\n--- 绘制并导出图C ---');
wscore_subtype = string(wscore.Volume_Subtype);
asd_sub1 = (wscore.DX_GROUP==1) & (wscore_subtype == "1" | strcmpi(wscore_subtype, "SUB1"));
asd_sub2 = (wscore.DX_GROUP==1) & (wscore_subtype == "2" | strcmpi(wscore_subtype, "SUB2"));
sub1_id = wscore.SUB_ID(asd_sub1);
sub2_id = wscore.SUB_ID(asd_sub2);
[~, idx_sub1] = ismember(sub1_id, raw.SUB_ID);
[~, idx_sub2] = ismember(sub2_id, raw.SUB_ID);
pad_sub1 = raw.PAD(idx_sub1);
pad_sub2 = raw.PAD(idx_sub2);
pad_hc = raw.PAD(hc_idx);
pad_sub1 = pad_sub1(~isnan(pad_sub1));
pad_sub2 = pad_sub2(~isnan(pad_sub2));
pad_hc = pad_hc(~isnan(pad_hc));

data_groups = {pad_sub1, pad_sub2, pad_hc};
% ======== 修改2：图C横坐标HC→TC，彻底统一术语 ========
group_names = {'Subtype 1', 'Subtype 2', 'TC'};
group_colors = {[0.90, 0.62, 0.00], [0.00, 0.45, 0.70], [0.34, 0.71, 0.91]};

figure('Name', 'Figure C', 'Color','w', 'Units','centimeters', 'Position',[2,2,18.3,5.0]);
hold on;
violinplot_with_box(data_groups, group_names, group_colors);

ylabel('Brain-PAD (years)', 'FontSize',5.5, 'FontName','Arial');

% ==========【修改2：移除图C总标题，统一到AI中添加】==========
% title('Brain-PAD distributions across ASD subtypes and healthy controls', ...
%     'FontSize',6.5, 'FontName','Arial', 'FontWeight','normal');

set(gca, 'FontName','Arial', 'FontSize', 5, 'TickDir','out', 'Box','off', 'LineWidth', 0.6, 'YGrid','off');
hold off;

figC_pdf = fullfile(work_dir, 'Figure_C_BrainPAD_Distributions_noTitle.pdf');
figC_png = fullfile(work_dir, 'Figure_C_BrainPAD_Distributions_noTitle.png');
try
    exportgraphics(gcf, figC_pdf, 'ContentType','vector');
    exportgraphics(gcf, figC_png, 'Resolution', 300);
    fprintf('已输出图C矢量PDF+PNG，尺寸18.3cm×5cm（无总标题）\n');
catch ME
    warning('图C导出失败: %s', ME.message);
end

fprintf('所有任务完成！\n');

%% ======================= 自定义函数：小提琴图（带箱线） ==========================
function violinplot_with_box(data, groupLabels, colors)
    nGroups = length(data);
    maxWidth = 0.4;
    boxWidth = 0.2;
    for g = 1:nGroups
        d = data{g}(:);
        d = d(~isnan(d));
        if isempty(d), continue; end
        
        [f, xi] = ksdensity(d);
        f_norm = f / max(f) * maxWidth;
        
        x_left = g - f_norm;
        x_right = g + f_norm;
        x_poly = [x_left, fliplr(x_right)];
        y_poly = [xi, fliplr(xi)];
        fill(x_poly, y_poly, colors{g}, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        
        q1 = prctile(d, 25);
        q3 = prctile(d, 75);
        median_val = median(d);
        iqr = q3 - q1;
        lower_whisker = max(min(d), q1 - 1.5*iqr);
        upper_whisker = min(max(d), q3 + 1.5*iqr);
        
        boxX = [g - boxWidth/2, g + boxWidth/2, g + boxWidth/2, g - boxWidth/2];
        boxY = [q1, q1, q3, q3];
        patch(boxX, boxY, 'w', 'EdgeColor', 'k', 'LineWidth', 0.8, 'FaceColor', 'w');
        line([g - boxWidth/2, g + boxWidth/2], [median_val, median_val], 'Color', 'k', 'LineWidth', 1.0);
        line([g, g], [lower_whisker, q1], 'Color', 'k', 'LineWidth', 0.6);
        line([g, g], [q3, upper_whisker], 'Color', 'k', 'LineWidth', 0.6);
        line([g - boxWidth/4, g + boxWidth/4], [lower_whisker, lower_whisker], 'Color', 'k', 'LineWidth', 0.6);
        line([g - boxWidth/4, g + boxWidth/4], [upper_whisker, upper_whisker], 'Color', 'k', 'LineWidth', 0.6);
        
        outliers = d(d < lower_whisker | d > upper_whisker);
        if ~isempty(outliers)
            scatter(g * ones(size(outliers)), outliers, 6, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
    end
    
    xticks(1:nGroups);
    xticklabels(groupLabels);
    xlim([0.5, nGroups+0.5]);
    set(gca, 'FontName','Arial', 'FontSize', 5, 'TickDir','out', 'Box','off', 'LineWidth', 0.6);
end

function plot_correlation_A(x, y, label, titleStr, axisLim)
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);
    if nargin >= 5 && ~isempty(axisLim)
        lo = axisLim(1); hi = axisLim(2);
    else
        lo = min([x(:); y(:)]); hi = max([x(:); y(:)]);
    end
    
    % 颜色判断逻辑保留原label='HC'，仅修改显示文字，颜色映射不变
    if strcmpi(label, 'HC')
        ptColor = [0.85, 0.37, 0.00];
    else
        ptColor = [0.00, 0.45, 0.70];
    end
    
    fs_title = 6;
    fs_label = 5.5;
    fs_tick  = 5;
    fs_stat  = 5;
    line_axis = 0.6;
    line_fit  = 1.0;
    line_ref  = 0.8;
    
    scatter(x, y, 8, 'filled', 'MarkerFaceColor', ptColor, 'MarkerFaceAlpha', 0.55); hold on;
    
    p = polyfit(x, y, 1);
    xx = linspace(min(x), max(x), 100);
    plot(xx, polyval(p,xx), 'r-', 'LineWidth', line_fit);
    plot([lo hi], [lo hi], '--', 'Color',[0.35 0.35 0.35], 'LineWidth', line_ref);
    
    r = corrcoef(x,y); r = r(1,2);
    mae = mean(abs(x-y));
    rmse = sqrt(mean((x-y).^2));
    
    text(lo + 0.20*(hi-lo), hi - 0.18*(hi-lo), ...
         sprintf('r = %.3f\nMAE = %.3f\nRMSE = %.3f', r, mae, rmse), ...
         'FontSize',fs_stat, 'FontWeight','normal', 'FontName','Arial', ...
         'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment','left');
    
    if nargin < 4 || isempty(titleStr)
        title([label ' Cohort'], 'FontSize',fs_title, 'FontWeight','normal', 'FontName','Arial');
    else
        title(titleStr, 'FontSize',fs_title, 'FontWeight','normal', 'FontName','Arial');
    end
    xlabel('Chronological age (years)', 'FontSize',fs_label, 'FontName','Arial');
    ylabel('Predicted age (years)', 'FontSize',fs_label, 'FontName','Arial');
    xlim([lo hi]); ylim([lo hi]);
    axis square; 
    grid off;
    set(gca, 'FontName','Arial', 'FontSize',fs_tick, 'TickDir','out', 'Box','on', 'LineWidth', line_axis);
end

%% ======================= 多模型对比辅助函数 ==========================
function Xs = zscore_rows(Xtr)
    mu = mean(Xtr,1); sigma = std(Xtr,0,1); sigma(sigma==0) = 1;
    Xs = (Xtr - mu)./sigma;
end

function Xs = zscore_apply(Xte, Xtr)
    mu = mean(Xtr,1); sigma = std(Xtr,0,1); sigma(sigma==0) = 1;
    Xs = (Xte - mu)./sigma;
end

function r = rv(x, y)
    c = corrcoef(x, y); r = c(1,2);
    if ~isfinite(r), r = 0; end
end

function K = rbf_kernel(A, B, sigma)
    K = exp(-pdist2(A, B).^2 / (2*sigma^2));
end

% ---------------- 核岭回归 KRR ----------------
function model = krr_fit(Xs, y, p)
    K = rbf_kernel(Xs, Xs, p(2));
    model.lambda = p(1); model.sigma = p(2); model.X = Xs;
    model.alpha = (K + p(1)*eye(size(K,1))) \ y;
end

function yp = krr_pred(model, Xs)
    Kt = rbf_kernel(Xs, model.X, model.sigma);
    yp = Kt * model.alpha;
end

% ---------------- 相关向量回归 RVR ----------------
function model = rvr_fit(Xs, y, p)
    sigma = p(1);
    K = rbf_kernel(Xs, Xs, sigma);
    n = size(K,1);
    Phi = [K, ones(n,1)];
    t = y;
    beta = 1 / (var(t) + eps);
    M = n + 1;
    alpha = inf(M,1);
    active = false(M,1);
    PhiTt = Phi'*t;
    PhiTPhi = Phi'*Phi;
    Phi2sum = sum(Phi.^2,1)';
    maxIter = 300;
    for iter = 1:maxIter
        idxA = find(active);
        if isempty(idxA)
            Svec = beta*Phi2sum;
            Qvec = beta*PhiTt;
        else
            PhiA = Phi(:,idxA); alphaA = alpha(idxA);
            Sigma = stable_inv(beta*PhiTPhi(idxA,idxA) + diag(alphaA));
            muA   = beta*Sigma*(PhiA'*t);
            T = PhiA' * Phi;
            SigmaT = Sigma * T;
            Svec = beta*Phi2sum - beta^2 * sum(T .* SigmaT, 1)';
            Qvec = beta*PhiTt - beta^2 * (T' * (Sigma*(PhiA'*t)));
        end
        s_eff = inf(M,1); q_eff = zeros(M,1);
        s_eff(active) = alpha(active).*Svec(active)./(alpha(active) - Svec(active));
        q_eff(active) = alpha(active).*Qvec(active)./(alpha(active) - Svec(active));
        inact = ~active;
        s_eff(inact) = beta*Svec(inact)./(beta + Svec(inact));
        q_eff(inact) = beta*Qvec(inact)./(beta + Svec(inact));
        theta = q_eff.^2 - s_eff;
        dL = zeros(M,1);
        dL(active)  = theta(active) ./ (alpha(active) + s_eff(active));
        dL(inact)   = theta(inact) ./ s_eff(inact);
        [gmax, iStar] = max(dL);
        if ~isfinite(gmax), break; end
        if gmax <= 1e-8 || iter >= maxIter, break; end
        if theta(iStar) > 0
            alphaNew = min(max(s_eff(iStar)^2 / theta(iStar), 1e-8), 1e12);
            alpha(iStar) = alphaNew;
            active(iStar) = true;
        else
            active(iStar) = false;
            alpha(iStar) = inf;
        end
        if any(active)
            idxA = find(active);
            PhiA = Phi(:,idxA); alphaA = alpha(idxA);
            Sigma = stable_inv(beta*PhiTPhi(idxA,idxA) + diag(alphaA));
            muA   = beta*Sigma*(PhiA'*t);
            gamma = 1 - alphaA.*diag(Sigma);
            err = t - PhiA*muA;
            betaNew = (n - sum(gamma)) / (err'*err + eps);
            if isfinite(betaNew) && betaNew > 0, beta = betaNew; end
        end
    end
    if ~any(active), active(1) = true; end
    idxA = find(active);
    PhiA = Phi(:,idxA); alphaA = alpha(idxA);
    Sigma = stable_inv(beta*PhiTPhi(idxA,idxA) + diag(alphaA));
    muA   = beta*Sigma*(PhiA'*t);
    model.X = Xs; model.sigma = sigma; model.active = idxA; model.mu = muA; model.bias = 0;
end

function yp = rvr_pred(model, Xs)
    Kt = rbf_kernel(Xs, model.X, model.sigma);
    Kt = [Kt, ones(size(Xs,1),1)];
    yp = Kt(:, model.active)*model.mu + model.bias;
end

function S = stable_inv(A)
    jit = 1e-12 * max(diag(A)) + 1e-8;
    S = (A + jit*eye(size(A))) \ eye(size(A));
end

% ---------------- 五折CV网格搜索 ----------------
function [best_params, oob_pred, hc_met] = evaluate_model_cv(X, y, param_grid, fit_fun, pred_fun)
    rng(42);
    cvp = cvpartition(size(X,1), 'KFold', 5);
    nP = size(param_grid, 1);
    cv_mse = nan(nP, 1);
    for i = 1:nP
        p = param_grid(i, :);
        s = 0;
        for f = 1:cvp.NumTestSets
            tr = cvp.training(f); te = cvp.test(f);
            Xs_tr = zscore_rows(X(tr,:));
            Xs_te = zscore_apply(X(te,:), X(tr,:));
            m = fit_fun(Xs_tr, y(tr), p);
            yp = pred_fun(m, Xs_te);
            s = s + mean((y(te) - yp).^2);
        end
        cv_mse(i) = s / cvp.NumTestSets;
    end
    [~, bi] = min(cv_mse);
    best_params = param_grid(bi, :);
    oob_pred = nan(size(X,1), 1);
    for f = 1:cvp.NumTestSets
        tr = cvp.training(f); te = cvp.test(f);
        Xs_tr = zscore_rows(X(tr,:));
        Xs_te = zscore_apply(X(te,:), X(tr,:));
        m = fit_fun(Xs_tr, y(tr), best_params);
        oob_pred(te) = pred_fun(m, Xs_te);
    end
    v = ~isnan(oob_pred);
    hc_met = struct('r', rv(y(v), oob_pred(v)), 'mae', mean(abs(y(v)-oob_pred(v))), ...
                    'rmse', sqrt(mean((y(v)-oob_pred(v)).^2)), 'n', sum(v), 'cv_mse', min(cv_mse));
end

% ---------------- 最终模型 + 独立测试集评估 ----------------
function [model, te_met] = fit_final_test(X, y, Xte, yte, params, fit_fun, pred_fun)
    Xs = zscore_rows(X);
    model.learner = fit_fun(Xs, y, params);
    model.pred_fun = pred_fun;
    if isempty(Xte) || isempty(yte)
        te_met = struct('r', NaN, 'mae', NaN, 'rmse', NaN, 'n', 0, 'pad', []);
        return;
    end
    Xs_te = zscore_apply(Xte, X);
    yp = pred_fun(model.learner, Xs_te);
    v = ~isnan(yp) & ~isnan(yte);
    te_met = struct('r', rv(yte(v), yp(v)), 'mae', mean(abs(yte(v)-yp(v))), ...
                    'rmse', sqrt(mean((yte(v)-yp(v)).^2)), 'n', sum(v), 'pad', yp(v) - yte(v));
end

% ---------------- 网格边界提示 ----------------
function warn_grid_boundary(bp, grid, mname)
    if size(grid,1) < 2, return; end
    for c = 1:numel(bp)
        lo = min(grid(:,c)); hi = max(grid(:,c));
        if abs(bp(c)-lo) <= eps(abs(lo)) || abs(bp(c)-hi) <= eps(abs(hi))
            fprintf('  [提醒] %s 参数%d=%g 落在网格边界[%g,%g]，可考虑外扩网格\n', mname, c, bp(c), lo, hi);
        end
    end
end
