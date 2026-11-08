%% ======================================================================
%  Fig4_05_ndm_permutation_test.m
%  NDM 拟合优度的非参数置换检验（补齐 Methods 2.7 所述、原实现缺失的环节）
%  ----------------------------------------------------------------------
%  背景（2026-09-18 审计）：原 Fig4_01 只输出 GOF 与 bootstrap 稳定性，输出
%  文件名即 NDM_Epicenter_Results_NoPermutation.xlsx，归档内没有任何置换检验
%  代码，而论文 Methods 2.7 声明"permutation tests were performed 5000 times"、
%  Results 3.3 声明 p < 0.0001。
%
%  本脚本：以球面旋转置换（保留空间自相关）扰动亚型平均异常图，重建 5000 次
%  零分布，检验最优拟合 Spearman R 的显著性：
%        p = (1 + #{R_null ≥ R_obs}) / (1 + n_perm)
%  注意：5000 次置换的最小可达 p = 1/5001 ≈ 2.0e-4，因此不可能出现 "p < 0.0001"。
%
%  性能优化：扩散预测模式 predicted(t, seed) 只依赖 MSN 与 t、与 atrophy 向量无关，
%  故对每个 t 预计算 68 个种子的预测模式秩向量（中心化+归一化），置换时仅需点积。
%            数学上与逐次调用 fit_ndm_epicenter 完全等价（脚本内含一致性断言）。
%
%  参数与 Fig4_01 严格一致：beta = 0.5、t ∈ linspace(0.1, 10, 50)、对称归一化拉普拉斯。
%  输入：ABIDE_All_Metrics_Subtypes_Summary_Wscore.xlsx、coord_lh_ct.csv、coord_rh_ct.csv
%  输出：NDM_Permutation_Results.xlsx
% ======================================================================
clear; clc;
tic;
SEED = 2024;
n_perm_target = 5000;
n_roi = 68;
beta = 0.5;                          % 与 Fig4_01 cfg.beta 一致
t_range = linspace(0.1, 10.0, 50);   % 与 Fig4_01 cfg.t_min/t_max/n_time_points 一致

%% ---- 0. 解析输入文件位置（允许从 code/ 目录直接运行）----
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cands = {here, fullfile(here,'..','data'), ...
         fullfile(here,'..','..','Figure1_Subtype_Identification','code'), ...
         fullfile(here,'..','..','Supplementary_Tools'), ...
         fullfile(here,'..','..','..')};
for i = 1:numel(cands), if exist(cands{i},'dir'), addpath(cands{i}); end, end  % 使 rotate_parcellation 等工具可用
resolve = @(fn) find_file(fn, cands);
FILE_W  = resolve('ABIDE_All_Metrics_Subtypes_Summary_Wscore.xlsx');
FILE_CL = resolve('coord_lh_ct.csv');
FILE_CR = resolve('coord_rh_ct.csv');
fprintf('输入：\n  %s\n  %s\n  %s\n', FILE_W, FILE_CL, FILE_CR);

%% ---- 1. MSN（严格复刻 Fig4_01）----
raw = readtable(FILE_W, 'VariableNamingRule', 'preserve');
raw.SUB_ID = string(raw.SUB_ID);
nv = raw.Properties.VariableNames(varfun(@isnumeric, raw, 'OutputFormat','uniform'));
for i = 1:numel(nv)
    v = nv{i};
    if ~strcmp(v, 'SUB_ID'), raw.(v)(raw.(v) <= -900) = NaN; end
end
metrics_all = {'Thickness','Curve','Abscurve','Area','Volume'};
tc = raw(raw.DX_GROUP == 2, :);
nf = nan(height(tc), n_roi, numel(metrics_all));
for m = 1:numel(metrics_all)
    cols = cellfun(@(x) sprintf('%s_%02d', metrics_all{m}, x), num2cell(1:n_roi), 'UniformOutput', false);
    nf(:,:,m) = tc{:, cols};
end
norm_f = zeros(size(nf));
for i = 1:size(nf,1)
    for r = 1:n_roi
        vec = squeeze(nf(i,r,:));
        if std(vec) > 0, norm_f(i,r,:) = (vec - mean(vec)) / std(vec); else, norm_f(i,r,:) = 0; end
    end
end
sum_MSN = zeros(n_roi); valid = 0;
for i = 1:size(norm_f,1)
    sm = corr(squeeze(norm_f(i,:,:))', 'type', 'Pearson');
    if ~any(isnan(sm(:))), sum_MSN = sum_MSN + sm; valid = valid + 1; end
end
conn = sum_MSN / valid;
conn(logical(eye(n_roi))) = 0;
conn = (conn + conn') / 2;
conn(conn < 0) = 0;
conn = conn / max(conn(:));
Dm = diag(1 ./ sqrt(sum(conn,2) + 1e-10));
L = eye(n_roi) - Dm * conn * Dm;      % 对称归一化拉普拉斯
fprintf('MSN 完成（TDC n=%d，有效 %d）\n', height(tc), valid);

%% ---- 2. 预计算各扩散时间点的预测模式秩矩阵 ----
[V, D] = eig((L + L')/2); lam = real(diag(D));
nT = numel(t_range);
A_rank = cell(nT, 1);
for it = 1:nT
    P = (V .* exp(-beta * lam * t_range(it))') * V';   % 第 i 列 = 以 ROI i 为种子的预测模式
    P = P ./ sum(P, 1);                                % 正值归一，不改变秩序
    R = zeros(size(P));
    for c = 1:size(P,2), R(:,c) = tiedrank_vec(P(:,c)); end
    Rc = R - mean(R, 1);
    A_rank{it} = Rc ./ vecnorm(Rc, 2, 1);
end
fprintf('预测模式预计算完成（%d 个时间点），用时 %.0f s\n', nT, toc);

%% ---- 3. 异常图与球面置换索引 ----
Ws = readtable(FILE_W);
Ws.SUB_ID = string(Ws.SUB_ID);
vol = cellfun(@(x) sprintf('Volume_%02d', x), num2cell(1:n_roi), 'UniformOutput', false);
asd_sel = (Ws.DX_GROUP == 1);
asdW = table2array(Ws(asd_sel, vol));
asd_lab = double(Ws.Volume_Subtype(asd_sel));
assert(numel(asd_lab) == size(asdW,1), 'mask 与数据行数不一致');

cl = readmatrix(FILE_CL); cr = readmatrix(FILE_CR);
cl = cl ./ vecnorm(cl, 2, 2); cr = cr ./ vecnorm(cr, 2, 2);
rng(SEED, 'twister');
perm_id = rotate_parcellation(cl, cr);
fprintf('置换索引 %s（使用前 %d 次，种子 %d）\n', mat2str(size(perm_id)), n_perm_target, SEED);

%% ---- 4. 观测值 + 零分布 ----
res_summary = cell(0, 8);
for s = 1:2
    mask = (asd_lab == s);
    av = mean(abs(asdW(mask,:)), 1, 'omitnan')';
    av = (av - min(av)) / (max(av) - min(av) + eps);

    [seed_o, t_o, r_obs, gof_all] = fit_fast(A_rank, t_range, av);
    [~, ord] = sort(gof_all, 'descend');
    % 一致性校验：与逐次特征分解的原始实现比对
    [s_ref, t_ref, r_ref, gof_ref] = fit_ref(L, av, t_range, beta);
    fprintf('\n[亚型 %d] n=%d | 观测: 震源 ROI %d | t=%.2f | R=%.4f | 前3名 %s\n', ...
        s, sum(mask), seed_o, t_o, r_obs, mat2str(ord(1:3)'));
    fprintf('      校验: 原始实现 ROI %d | t=%.2f | R=%.4f → max|ΔR|=%.2e, max|ΔGOF|=%.2e\n', ...
        s_ref, t_ref, r_ref, abs(r_obs - r_ref), max(abs(gof_all - gof_ref)));
    assert(abs(r_obs - r_ref) < 1e-9 && max(abs(gof_all - gof_ref)) < 1e-9, ...
        '快速实现与原始实现不一致，请检查 A_rank 预计算');

    r_null = zeros(n_perm_target, 1);
    for b = 1:n_perm_target
        [~, ~, r_null(b)] = fit_fast(A_rank, t_range, av(perm_id(:,b)));
    end
    p_perm = (1 + sum(r_null >= r_obs)) / (1 + n_perm_target);
    fprintf('   ★ 置换检验: p = %.5f（零分布中位 %.3f，最大 %.3f）\n', p_perm, median(r_null), max(r_null));

    res_summary = [res_summary; {s, sum(mask), seed_o, t_o, r_obs, p_perm, median(r_null), max(r_null)}]; %#ok<AGROW>
    writetable(cell2table(res_summary, 'VariableNames', ...
        {'Subtype','N','Best_ROI','Best_t','Observed_R','P_perm','Null_median','Null_max'}), ...
        'NDM_Permutation_Results.xlsx');
end
disp(cell2table(res_summary, 'VariableNames', ...
    {'Subtype','N','Best_ROI','Best_t','Observed_R','P_perm','Null_median','Null_max'}));
fprintf('完成：NDM_Permutation_Results.xlsx（总用时 %.0f s）\n', toc);

%% ---- 局部函数 ----
function [best_seed, best_t, best_corr, gof_all] = fit_fast(A_rank, t_range, atrophy_vec)
    % 快速实现：使用预计算的预测模式秩矩阵
    br = tiedrank_vec(atrophy_vec(:)); br = br - mean(br);
    bn = br / norm(br);
    nT = numel(t_range); n_roi = numel(atrophy_vec);
    M = zeros(nT, n_roi);
    for it = 1:nT
        M(it, :) = (A_rank{it}' * bn)';
    end
    gof_all = max(M, [], 1)';
    [best_corr, idx] = max(M(:));
    [it_best, best_seed] = ind2sub(size(M), idx);
    best_t = t_range(it_best);
end

function [best_seed, best_t, best_corr, gof_all] = fit_ref(L, atrophy_vec, t_range, beta)
    % 原始逐次实现（与 Fig4_01 完全一致），仅用于校验快速实现
    n_roi = size(L,1);
    L = (L + L')/2;
    [V, D] = eig(L); lam = real(diag(D));
    b = tiedrank_vec(atrophy_vec(:)); b = b - mean(b); nb = norm(b);
    best_corr = -Inf; best_seed = -1; best_t = -1; gof_all = zeros(n_roi,1);
    for seed_idx = 1:n_roi
        Vt_seed = V(seed_idx,:)';
        max_seed_r = -Inf;
        for t = t_range
            predicted = V * (exp(-beta * lam * t) .* Vt_seed);
            predicted = predicted / sum(predicted);
            a = tiedrank_vec(predicted); a = a - mean(a);
            r = (a' * b) / (norm(a) * nb);
            if r > max_seed_r, max_seed_r = r; end
            if r > best_corr, best_corr = r; best_seed = seed_idx; best_t = t; end
        end
        gof_all(seed_idx) = max_seed_r;
    end
end

function r = tiedrank_vec(x)
    n = numel(x); [xs, ord] = sort(x); r = zeros(n,1); r(ord) = 1:n;
    i = 1;
    while i <= n
        j = i;
        while j < n && xs(j+1) == xs(i), j = j + 1; end
        if j > i, r(ord(i:j)) = mean(i:j); end
        i = j + 1;
    end
end

function p = find_file(fn, cands)
    % 在候选目录中查找输入文件；找不到则报错并列出已尝试路径
    p = '';
    for i = 1:numel(cands)
        c = fullfile(cands{i}, fn);
        if exist(c, 'file'), p = c; return; end
    end
    error('未找到 %s；已尝试：\n%s', fn, strjoin(cands, '\n'));
end
