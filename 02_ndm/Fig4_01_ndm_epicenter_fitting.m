%% ABIDE 多参数形态学相似性网络(MSN)与网络扩散模型(NDM)震源分析系统
% 【修改】删除了置换检验，添加了前7名拟合优度展示
clear; clc; close all;
warning('off', 'all');

fprintf('============================================================\n');
fprintf('   🧠 基于 MSN 网络扩散模型 (NDM) 震源分析系统 (精简版，无置换检验)\n');
fprintf('============================================================\n');

% =========================================================================
% --- 1. 参数配置 ---
% =========================================================================
cfg.n_roi = 68;                  % ROI数量 (DK-68图谱)
cfg.beta = 0.5;                  % 扩散速率常数
cfg.t_min = 0.1;                 % 最小传播时间
cfg.t_max = 10.0;                % 最大传播时间
cfg.n_time_points = 50;          % 时间步长离散点数

% Bootstrap 开关（保留）
cfg.run_bootstrap = true;       
cfg.n_bootstrap = 1000;         % 按要求改为 1000 次重采样（原为 100）

cfg.draw_figures = false;        % 关闭绘图

cfg.metrics = {'Volume'};        % 只分析 Volume

% 震源集群筛选配置
cfg.top_mode = 'quantile';        % 'quantile' 或 'topN'
cfg.top_thresh = 0.85;            % 分位数阈值
cfg.topN = 12;                    % 固定前N个

% ===================== 路径设置 =====================
% ---- Repository paths (no edit needed if you keep this folder layout) ----
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(script_dir);            % repository root
data_file    = fullfile(repo_root, 'data', 'ABIDE_All_Metrics_Subtypes_Summaryx.xlsx');
save_gof_dir = fullfile(repo_root, 'output', 'NDM_GOF_Maps');
if ~exist(save_gof_dir, 'dir'), mkdir(save_gof_dir); end

if ~exist(data_file, 'file')
    error('【文件缺失】: 路径下找不到数据集文件 → %s', data_file);
end
fprintf('✅ 文件加载成功：%s\n', data_file);

% =========================================================================
% --- 2. 构建 MSN 网络（完全保留）---
% =========================================================================
fprintf('\n[1/5] 正在从健康对照组 (TC, DX_GROUP==2) 提取多参数构建 MSN... \n');
dataT = readtable(data_file, 'VariableNamingRule', 'preserve');

tc_mask = (dataT.DX_GROUP == 2);
if sum(tc_mask) == 0, error('未检测到 TC 组！'); end
T_tc = dataT(tc_mask, :);
n_tc = height(T_tc);
fprintf('      健康对照样本量: %d 例\n', n_tc);

% 提取5种形态指标
full_metrics_all = {'Volume','Area','Thickness','Curve','Abscurve'};
raw_features = zeros(n_tc, cfg.n_roi, length(full_metrics_all));
for m = 1:length(full_metrics_all)
    target = full_metrics_all{m};
    roi_cols = cellfun(@(x) sprintf('%s_%02d', target, x), num2cell(1:cfg.n_roi), 'UniformOutput', false);
    raw_features(:, :, m) = T_tc{:, roi_cols};
end

% 跨特征标准化
fprintf('      ⚙️ 特征间标准化... \n');
norm_features = zeros(n_tc, cfg.n_roi, length(full_metrics_all));
for i = 1:n_tc
    for r = 1:cfg.n_roi
        vec = squeeze(raw_features(i, r, :));
        if std(vec) > 0
            norm_features(i, r, :) = (vec - mean(vec)) / std(vec);
        else
            norm_features(i, r, :) = 0;
        end
    end
end

% 组平均MSN
sum_MSN = zeros(cfg.n_roi, cfg.n_roi);
valid_sub = 0;
for i = 1:n_tc
    sub_matrix = squeeze(norm_features(i, :, :));
    sub_msn = corr(sub_matrix', 'type', 'Pearson');
    if ~any(isnan(sub_msn(:)))
        sum_MSN = sum_MSN + sub_msn;
        valid_sub = valid_sub + 1;
    end
end
conn_matrix = sum_MSN / valid_sub;
conn_matrix(logical(eye(cfg.n_roi))) = 0;
conn_matrix = (conn_matrix + conn_matrix') / 2;
conn_matrix(conn_matrix < 0) = 0;
conn_matrix = conn_matrix / max(conn_matrix(:));
fprintf('      ✅ MSN 网络构建完成 (%dx%d)\n', cfg.n_roi, cfg.n_roi);

% 拉普拉斯矩阵
degrees = sum(conn_matrix, 2);
D_inv_sqrt = diag(1.0 ./ sqrt(degrees + 1e-10));
A_norm = D_inv_sqrt * conn_matrix * D_inv_sqrt;
L_matrix = eye(cfg.n_roi) - A_norm;

% =========================================================================
% --- 3. 对每个亚型进行分析 ---
% =========================================================================
t_range = linspace(cfg.t_min, cfg.t_max, cfg.n_time_points);
roi_labels = cellfun(@(x) sprintf('ROI_%02d', x), num2cell(1:cfg.n_roi), 'UniformOutput', false);
res_summary = cell(0, 8);   % 结果表，少了 p 值列

for m = 1:length(cfg.metrics)
    target = cfg.metrics{m};
    fprintf('\n------------------------------------------------------------\n');
    fprintf('🚀 分析指标: [%s]\n', target);
    fprintf('------------------------------------------------------------\n');
    
    roi_cols = cellfun(@(x) sprintf('%s_%02d', target, x), num2cell(1:cfg.n_roi), 'UniformOutput', false);
    roi_data = dataT{:, roi_cols};
    
    subtype_label = [target, '_Subtype'];
    subtypes = dataT.(subtype_label);
    if iscell(subtypes) || isstring(subtypes) || iscategorical(subtypes)
        subtypes_num = zeros(length(subtypes), 1);
        subtypes_num(contains(string(subtypes), '1')) = 1;
        subtypes_num(contains(string(subtypes), '2')) = 2;
        subtypes = subtypes_num;
    end
    
    unique_sub = [1,2];
    
    for s = 1:length(unique_sub)
        s_val = unique_sub(s);
        pts_mask = (dataT.DX_GROUP == 1 & subtypes == s_val);
        n_pts = sum(pts_mask);
        if n_pts < 5
            fprintf('  [亚型 %d] 样本量不足 (%d 例)，跳过。\n', s_val, n_pts);
            continue;
        end
        
        % 计算亚型平均异常幅度
        individual_abs = abs(roi_data(pts_mask, :));
        atrophy_vec = mean(individual_abs, 1, 'omitnan')';
        atrophy_vec = (atrophy_vec - min(atrophy_vec)) / (max(atrophy_vec) - min(atrophy_vec) + eps);
        
        % ---------- NDM拟合 ----------
        [best_seed, best_t, best_corr, ndm_gof_allroi] = fit_ndm_epicenter(L_matrix, atrophy_vec, t_range, cfg.beta);
        fprintf('  🎯 [Subtype %d] 全局最优震源: %d 号 | t=%.2f | R=%.3f\n', s_val, best_seed, best_t, best_corr);
        
        % ---------- 展示前7名拟合优度 ----------
        [sorted_gof, sorted_idx] = sort(ndm_gof_allroi, 'descend');
        top7_idx = sorted_idx(1:7);
        top7_gof = sorted_gof(1:7);
        fprintf('  📊 前7名拟合优度:\n');
        for k = 1:7
            fprintf('      第%d名: ROI %3d  (GOF = %.4f)\n', k, top7_idx(k), top7_gof(k));
        end
        
        % 保存前7名到Excel（可选）
        top7_table = table((1:7)', top7_idx, top7_gof, ...
            'VariableNames', {'Rank', 'ROI', 'GOF'});
        top7_excel = fullfile(save_gof_dir, sprintf('Top7_GOF_%s_Subtype%d.xlsx', target, s_val));
        writetable(top7_table, top7_excel);
        
        % ---------- 震源集群筛选 ----------
        [epicenter_cluster, cutoff_val] = extract_top_epicenters(ndm_gof_allroi, cfg);
        cluster_str = strjoin(cellstr(num2str(epicenter_cluster(:))), ',');
        fprintf('             ┗━ 震源集群 ROIs: [%s] | 阈值=%.3f\n', cluster_str, cutoff_val);
        
        % 保存全脑GOF
        gof_save_name = sprintf('GOF_Map_%s_Subtype%d.mat', target, s_val);
        gof_fullpath = fullfile(save_gof_dir, gof_save_name);
        save(gof_fullpath, 'ndm_gof_allroi', 'epicenter_cluster', 'cutoff_val');
        
        % 导出全脑GOF Excel
        roi_num_list = (1:cfg.n_roi)';
        full_gof_table = table(roi_num_list, ndm_gof_allroi, 'VariableNames', {'ROI编号', '拟合优度GOF'});
        excel_gof_name = sprintf('AllROI_FullGOF_%s_Subtype%d.xlsx', target, s_val);
        writetable(full_gof_table, fullfile(save_gof_dir, excel_gof_name));
        
        % ---------- Bootstrap稳定性 ----------
        stability_score = NaN;
        if cfg.run_bootstrap
            % 【2026-09-19 修正】原实现误在 68 个 ROI 上重采样并重排 MSN 矩阵；
            % 现按论文 2.7 口径改为【被试级】有放回重采样（MSN 固定不变）。
            [stability_score, ~] = bootstrap_stability(roi_data(pts_mask, :), L_matrix, ...
                                                       t_range, cfg.beta, cfg.n_bootstrap);
            fprintf('             ┗━ Bootstrap稳定性: %.2f%%\n', stability_score * 100);
        end
        
        % 汇总（无 p 值）
        res_summary = [res_summary; {target, s_val, n_pts, best_seed, cluster_str, best_corr, cutoff_val, stability_score}];
    end
end

% =========================================================================
% --- 4. 导出汇总报表 ---
% =========================================================================
fprintf('\n[5/5] 正在保存汇总结果...\n');
resTable = cell2table(res_summary, 'VariableNames', ...
    {'Metric', 'Subtype', 'SampleSize', 'Single_Best_Epicenter', 'Epicenter_Cluster_ROIs', ...
     'Max_Spearman_R', 'GOF_Cutoff_Threshold', 'Bootstrap_Stability'});
writetable(resTable, fullfile(repo_root, 'output', 'NDM_Epicenter_Results_NoPermutation.xlsx'));
fprintf('      Summary table saved to output/NDM_Epicenter_Results_NoPermutation.xlsx\n');

%% 加载并查看MSN矩阵（可选）
fprintf('\n=============================================\n');
fprintf('📌 MSN 矩阵摘要\n');
fprintf('=============================================\n');
fprintf('MSN矩阵尺寸：%d × %d\n', size(conn_matrix,1), size(conn_matrix,2));
fprintf('MSN矩阵最大值：%.4f，最小值：%.4f\n', max(conn_matrix(:)), min(conn_matrix(:)));

% 可选：绘制热力图（如需要可取消注释）
% figure; imagesc(conn_matrix); colorbar; axis square; title('MSN矩阵');

fprintf('\n==================== 🎉 分析完成 ====================\n');

% =========================================================================
% --- 辅助函数定义 ---
% =========================================================================

function [best_seed, best_t, best_corr, ndm_gof_allroi] = fit_ndm_epicenter(L, atrophy_vec, t_range, beta)
    n_roi = size(L,1);
    L = (L + L') / 2;                 % 数值对称化，保证 eig 稳定（L 理论上对称）
    [V, D] = eig(L);                  % 对称矩阵特征分解：e^{-beta*L*t} = V*diag(e^{-beta*lam*t})*V'
    lam = real(diag(D));
    % 预先计算目标向量的秩（Spearman 秩相关 = 秩上的 Pearson 相关）
    b = rank_vec(atrophy_vec(:));
    b = b - mean(b); nb = norm(b);
    best_corr = -Inf; best_seed = -1; best_t = -1;
    ndm_gof_allroi = zeros(n_roi,1);
    for seed_idx = 1:n_roi
        Vt_seed = V(seed_idx,:)';     % V' * e_seed（单位种子向量）
        max_seed_r = -Inf;
        for t = t_range
            % 等价于 expm(-beta*L*t)*seed_vec，但为 O(n^2)
            predicted = V * (exp(-beta * lam * t) .* Vt_seed);
            predicted = predicted / sum(predicted);
            a = rank_vec(predicted);
            a = a - mean(a);
            r = (a' * b) / (norm(a) * nb);
            if r > max_seed_r
                max_seed_r = r;
            end
            if r > best_corr
                best_corr = r;
                best_seed = seed_idx;
                best_t = t;
            end
        end
        ndm_gof_allroi(seed_idx) = max_seed_r;
    end
end

function r_out = rank_vec(x)
    % 快速秩变换（平均秩处理并列值），替代 tiedrank/corr 的高开销调用
    n = numel(x);
    [xs, ord] = sort(x);
    r_out = zeros(n,1);
    r_out(ord) = 1:n;
    % 处理并列值：同一数值取平均秩
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i)
            j = j + 1;
        end
        if j > i
            r_out(ord(i:j)) = mean(i:j);
        end
        i = j + 1;
    end
end

function [epicenter_roi_list, cutoff] = extract_top_epicenters(gof_vec, cfg)
    switch cfg.top_mode
        case 'quantile'
            cutoff = quantile(gof_vec, cfg.top_thresh);
            epicenter_roi_list = find(gof_vec >= cutoff);
        case 'topN'
            [~, idx_sort] = sort(gof_vec, 'descend');
            epicenter_roi_list = idx_sort(1:cfg.topN);
            cutoff = gof_vec(idx_sort(cfg.topN));
    end
    epicenter_roi_list = sort(epicenter_roi_list);
end

function [stability_score, seed_freqs] = bootstrap_stability(roi_mat, L_matrix, t_range, beta, n_boot)
    % 【2026-09-19 修正】被试级 Bootstrap（与论文 2.7 描述一致）：
    %   每次迭代有放回抽取【被试】→ 重算亚型平均 |W| 异常图 → MSN/拉普拉斯固定不变
    %   → 重跑网络扩散拟合，统计各脑区被选为"最优震源"的频率。
    %   原实现错误地对 68 个 ROI 重采样并同时重排了 MSN 矩阵，无统计学含义。
    % 加速：扩散预测模式只依赖 MSN 与扩散时间 t（与异常图无关），故预先算好
    %       每个 t、每个种子的预测模式秩向量，迭代时仅做点积；与逐次调用
    %       fit_ndm_epicenter 数学等价。
    rng(2024, 'twister');                 % 固定种子，保证可复现
    n_roi = size(L_matrix, 1);
    [V, D] = eig((L_matrix + L_matrix') / 2); lam = real(diag(D));
    nT = numel(t_range);
    A_rank = cell(nT, 1);
    for it = 1:nT
        P = (V .* exp(-beta * lam * t_range(it))') * V';
        P = P ./ sum(P, 1);
        R = zeros(n_roi);
        for c = 1:n_roi, R(:, c) = rank_vec(P(:, c)); end
        Rc = R - mean(R, 1);
        A_rank{it} = Rc ./ vecnorm(Rc, 2, 1);
    end
    Xabs = abs(roi_mat);
    n_sub = size(Xabs, 1);
    seed_counts = zeros(n_roi, 1);
    for b = 1:n_boot
        bi = randi(n_sub, [n_sub, 1]);                     % ← 被试级有放回重采样
        av = mean(Xabs(bi, :), 1, 'omitnan')';
        av = (av - min(av)) / (max(av) - min(av) + eps);    % 与主流程一致的 0-1 归一
        br = rank_vec(av); br = br - mean(br); br = br / norm(br);
        gof = -inf(n_roi, 1);
        for it = 1:nT
            gof = max(gof, A_rank{it}' * br);               % MSN 固定，不重采样
        end
        [~, bs] = max(gof);
        seed_counts(bs) = seed_counts(bs) + 1;
    end
    seed_freqs = seed_counts / n_boot;
    stability_score = max(seed_freqs);
end