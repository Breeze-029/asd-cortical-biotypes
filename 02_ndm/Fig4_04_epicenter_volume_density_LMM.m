%% =========================================================================
%  震源区 vs 非震源区：皮层体积绝对偏差 核密度验证（Nature 投稿配图风格）
%  输入：ABIDE_All_Metrics_Subtypes_Summary_Wscore.xlsx（Volume W-score）
%  主检验：线性混合模型 LMM（AbsDev ~ RegionType + (1|Subject)），避免伪重复
%  输出：Figure_Volume_Epicenter_density.pdf / .png / .tiff
%  规范：Arial；183 × 65 mm（双栏整宽，横向两面板并列）；
%        细线（0.5–0.75 pt）；矢量 PDF + 600 dpi 位图；配色与本文件夹其它配图统一
% =========================================================================
clear; clc; close all;
warning('off', 'all');

%% ========================= 1. 读取数据 =========================
% ---- Repository paths (no edit needed if you keep this folder layout) ----
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(script_dir);            % repository root
filePath = fullfile(repo_root, 'data', 'ABIDE_All_Metrics_Subtypes_Summary_Wscore.xlsx');
data = readtable(filePath, 'VariableNamingRule', 'preserve');

%% ========================= 2. 提取体积指标与分组 =========================
volume_cols = cellfun(@(x) sprintf('Volume_%02d', x), num2cell(1:68), 'UniformOutput', false);
volume_data = data{:, volume_cols};
subtype     = data.Volume_Subtype;    % 体积亚型
dx_group    = data.DX_GROUP;          % 1 = ASD, 2 = TDC

% 仅保留 ASD 患者
patient_mask = dx_group == 1;
volume_data  = volume_data(patient_mask, :);
subtype      = subtype(patient_mask);

%% ========================= 3. NDM 识别的震源区（ROI 编号）=========================
% 与正文一致的候选震源集群 = GOF 排名前 10 的脑区（对应 quantile 0.85、阈值 0.675/0.199）
% Subtype 1: 最优震源 ROI 18 (right medial orbitofrontal)，集群 lateral frontotemporal-insula
epi_roi_sub1 = [15, 16, 18, 19, 25, 26, 55, 56, 58, 59];
% Subtype 2: 最优震源 ROI 62 (right temporal pole)，集群 anterior temporal-limbic
epi_roi_sub2 = [14, 22, 24, 25, 26, 35, 43, 44, 61, 62];
all_roi = 1:68;

%% ========================= 4. 提取绝对体积偏差（|W-score|）=========================
data1 = volume_data(subtype == 1, :);
data2 = volume_data(subtype == 2, :);

epi_vals1 = [];  non_vals1 = [];
for i = 1:size(data1,1)
    row = abs(data1(i,:));
    epi_vals1 = [epi_vals1; row(epi_roi_sub1)'];
    non_vals1 = [non_vals1; row(setdiff(all_roi, epi_roi_sub1))'];
end

epi_vals2 = [];  non_vals2 = [];
for i = 1:size(data2,1)
    row = abs(data2(i,:));
    epi_vals2 = [epi_vals2; row(epi_roi_sub2)'];
    non_vals2 = [non_vals2; row(setdiff(all_roi, epi_roi_sub2))'];
end

%% ========================= 5. 主检验：线性混合模型（LMM）=========================
% 被试作为随机截距：同一被试贡献 7 个震源区值与 61 个非震源区值，观测非独立，
% 观测级 Wilcoxon（n=3171 vs 27633）会因伪重复使 p 值不可解释。
[p_val_sub1, beta_sub1, se_sub1, t_sub1] = lmm_epicenter_test(data1, epi_roi_sub1);
[p_val_sub2, beta_sub2, se_sub2, t_sub2] = lmm_epicenter_test(data2, epi_roi_sub2);

% 对照：观测级 Wilcoxon（仅作记录，不作主结论）
w_sub1 = ranksum(epi_vals1, non_vals1);
w_sub2 = ranksum(epi_vals2, non_vals2);

fprintf('==================== Epicenter volume validation ====================\n');
fprintf('Subtype 1: n = %d | median Epi = %.4f, Non = %.4f\n', ...
    size(data1,1), median(epi_vals1), median(non_vals1));
fprintf('  LMM: beta(Non-Epi) = %.4f (SE = %.4f), t = %.2f, p = %.4e   [main test]\n', ...
    beta_sub1, se_sub1, t_sub1, p_val_sub1);
fprintf('  Wilcoxon (observation-level, reference) p = %.4e | Cliff''s delta = %.4f\n', ...
    w_sub1, cliffs_delta(epi_vals1, non_vals1));
fprintf('Subtype 2: n = %d | median Epi = %.4f, Non = %.4f\n', ...
    size(data2,1), median(epi_vals2), median(non_vals2));
fprintf('  LMM: beta(Non-Epi) = %.4f (SE = %.4f), t = %.2f, p = %.4e   [main test]\n', ...
    beta_sub2, se_sub2, t_sub2, p_val_sub2);
fprintf('  Wilcoxon (observation-level, reference) p = %.4e | Cliff''s delta = %.4f\n', ...
    w_sub2, cliffs_delta(epi_vals2, non_vals2));
fprintf('=====================================================================\n');

%% ========================= 6. 绘图（Nature 风格，两面板同尺度）=========================
font_name = 'Arial';
fs_ax    = 7;      % 坐标轴标签
fs_tick  = 6;      % 刻度数字
fs_leg   = 6;      % 图例
fs_panel = 7;      % 面板标题
lw_curve = 0.75;   % 密度曲线线宽
lw_axis  = 0.5;    % 坐标轴/箱体线宽

% 配色：与 Fig2_C / Fig2_E 统一（暖红 = 震源区；冷蓝 = 非震源区）
c_epi = [0.78 0.45 0.45];
c_non = [0.36 0.54 0.66];
% 填充色：预先与白底按 30% 透明度混合
% （矢量 PDF 不支持 alpha 透明，预混色可保证 PDF 与屏幕外观完全一致）
alpha_fill = 0.30;
c_epi_fill = alpha_fill * c_epi + (1 - alpha_fill);
c_non_fill = alpha_fill * c_non + (1 - alpha_fill);

% 两个面板共用横轴范围，便于直接比较分布
x_all = [epi_vals1; non_vals1; epi_vals2; non_vals2];
x_all = x_all(~isnan(x_all));
xlim_shared = [0, ceil(max(x_all) * 10.2) / 10];

% 画布 = 投稿物理尺寸（双栏整宽 183 mm × 高 65 mm），两个分布图左右并列
fig = figure('Units','centimeters','Position',[1 1 18.3 6.5],'Color','w', ...
             'PaperUnits','centimeters','PaperSize',[18.3 6.5], ...
             'PaperPosition',[0 0 18.3 6.5],'PaperPositionMode','manual');

tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl);
plot_density_curves(ax1, epi_vals1, non_vals1, xlim_shared, 'Subtype 1', p_val_sub1, ...
    c_epi, c_non, c_epi_fill, c_non_fill, font_name, fs_ax, fs_tick, fs_leg, fs_panel, ...
    lw_curve, lw_axis);

ax2 = nexttile(tl);
plot_density_curves(ax2, epi_vals2, non_vals2, xlim_shared, 'Subtype 2', p_val_sub2, ...
    c_epi, c_non, c_epi_fill, c_non_fill, font_name, fs_ax, fs_tick, fs_leg, fs_panel, ...
    lw_curve, lw_axis);

%% ========================= 7. 导出（矢量 PDF + 600 dpi 位图）=========================
out_dir  = fullfile(repo_root, 'output');
pdf_file = fullfile(out_dir, 'Figure_Volume_Epicenter_density.pdf');
png_file = fullfile(out_dir, 'Figure_Volume_Epicenter_density.png');
tif_file = fullfile(out_dir, 'Figure_Volume_Epicenter_density.tiff');

% 用 print + PaperSize 精确控制页面物理尺寸（183 × 65 mm，横向两面板并列）
% 填充使用预混浅色，painters 矢量输出不涉及透明，PDF 与屏幕外观一致
print(fig, pdf_file, '-dpdf', '-vector');
print(fig, png_file, '-dpng', '-r600');
print(fig, tif_file, '-dtiff', '-r600');

fprintf('导出完成：%s（183 × 65 mm，横向两面板）\n', pdf_file);

%% ========================= 辅助函数 =========================
function plot_density_curves(ax, data_epi, data_non, xlims, panel_title, p_val, ...
    c_epi, c_non, c_epi_fill, c_non_fill, font_name, fs_ax, fs_tick, fs_leg, fs_panel, ...
    lw_curve, lw_axis)
% 单面板：两组核密度曲线 + 底部迷你箱线 + 显著性标记（Nature 风格）

    [f1, xi1] = ksdensity(data_epi, 'NumPoints', 300);
    [f2, xi2] = ksdensity(data_non, 'NumPoints', 300);
    y_max = max([f1(:); f2(:)]);

    hold(ax, 'on');

    % 核密度填充（预混浅色，保证矢量 PDF 外观与屏幕一致）+ 轮廓线
    fill(ax, [xi1, fliplr(xi1)], [f1, zeros(size(f1))], c_epi_fill, ...
        'FaceColor', c_epi_fill, 'EdgeColor', 'none');
    fill(ax, [xi2, fliplr(xi2)], [f2, zeros(size(f2))], c_non_fill, ...
        'FaceColor', c_non_fill, 'EdgeColor', 'none');
    plot(ax, xi1, f1, 'Color', c_epi, 'LineWidth', lw_curve);
    plot(ax, xi2, f2, 'Color', c_non, 'LineWidth', lw_curve);

    % 底部迷你箱线（四分位 + 中位数），位置随密度量级缩放
    q1 = quantile(data_epi, [0.25 0.50 0.75]);
    q2 = quantile(data_non, [0.25 0.50 0.75]);
    y_base = -0.10 * y_max;
    h_box  =  0.05 * y_max;
    rectangle(ax, 'Position', [q1(1), y_base, q1(3)-q1(1), h_box], ...
        'FaceColor', c_epi, 'EdgeColor', 'k', 'LineWidth', lw_axis);
    rectangle(ax, 'Position', [q2(1), y_base, q2(3)-q2(1), h_box], ...
        'FaceColor', c_non, 'EdgeColor', 'k', 'LineWidth', lw_axis);
    line(ax, [q1(2) q1(2)], [y_base, y_base + h_box], 'Color', 'k', 'LineWidth', lw_axis);
    line(ax, [q2(2) q2(2)], [y_base, y_base + h_box], 'Color', 'k', 'LineWidth', lw_axis);

    % 显著性星号（基于 LMM p）
    if     p_val < 0.001, sig = '***';
    elseif p_val < 0.01,  sig = '**';
    elseif p_val < 0.05,  sig = '*';
    else,                 sig = 'n.s.';
    end
    text(ax, mean([q1(2), q2(2)]), y_max * 1.04, sig, ...
        'FontName', font_name, 'FontSize', fs_panel + 1, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    % 坐标轴与标签
    xlim(ax, xlims);
    ylim(ax, [-0.17 * y_max, y_max * 1.16]);
    xlabel(ax, 'Absolute volume deviation', 'FontName', font_name, 'FontSize', fs_ax);
    ylabel(ax, 'Density', 'FontName', font_name, 'FontSize', fs_ax);

    % 面板标题（第二行标注 LMM 统计，便于正文直接引用）
    title(ax, {panel_title, sprintf('LMM p = %.1e', p_val)}, ...
        'FontName', font_name, 'FontSize', fs_panel, 'FontWeight', 'bold');

    legend(ax, {'Epicenter', 'Non-epicenter'}, 'Location', 'northeast', ...
        'FontName', font_name, 'FontSize', fs_leg, 'Box', 'off');

    % Nature 风格坐标轴：仅保留左/下轴线（box off），刻度朝外，细线
    set(ax, 'FontName', font_name, 'FontSize', fs_tick, 'TickDir', 'out', ...
        'Box', 'off', 'LineWidth', lw_axis, 'Layer', 'top');
    grid(ax, 'off');
end

function delta = cliffs_delta(x, y)
% Cliff's delta 效应量（观测级，仅作参考）
    x = x(~isnan(x));
    y = y(~isnan(y));
    n1 = numel(x);  n2 = numel(y);
    greater = 0;    smaller = 0;
    for i = 1:n1
        greater = greater + sum(x(i) > y);
        smaller = smaller + sum(x(i) < y);
    end
    delta = (greater - smaller) / (n1 * n2);
end

function [p_val, beta, se, tstat] = lmm_epicenter_test(data_sub, epi_rois)
% 线性混合模型：AbsDev ~ RegionType + (1|Subject)
% 被试作为随机截距，正确建模"同一被试多个 ROI 观测"的非独立性（避免伪重复）
    n_sub    = size(data_sub, 1);
    all_roi  = 1:size(data_sub, 2);
    non_rois = setdiff(all_roi, epi_rois);
    n_epi = numel(epi_rois);
    n_non = numel(non_rois);
    n = n_sub * (n_epi + n_non);
    subj    = zeros(n, 1);
    regtype = cell(n, 1);
    absdev  = zeros(n, 1);
    k = 0;
    for i = 1:n_sub
        rowv = abs(data_sub(i, :));
        for r = epi_rois
            k = k + 1; subj(k) = i; regtype{k} = 'Epicenter';     absdev(k) = rowv(r);
        end
        for r = non_rois
            k = k + 1; subj(k) = i; regtype{k} = 'Non-epicenter'; absdev(k) = rowv(r);
        end
    end
    T = table(subj, categorical(regtype), absdev, ...
              'VariableNames', {'Subject', 'RegionType', 'AbsDev'});
    lme   = fitlme(T, 'AbsDev ~ RegionType + (1|Subject)');
    idx   = find(contains(string(lme.Coefficients.Name), 'Non'));
    p_val = lme.Coefficients.pValue(idx);
    beta  = lme.Coefficients.Estimate(idx);
    se    = lme.Coefficients.SE(idx);
    tstat = lme.Coefficients.tStat(idx);
end
