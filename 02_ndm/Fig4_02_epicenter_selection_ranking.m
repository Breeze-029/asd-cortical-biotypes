%% ========================================================================
%  基于NDM拟合优度(GOF)筛选疾病震源(Epicenter) - 完整工具
%  输入：data/Fig4_allROI_GOF_volume.xlsx（三列：ROI, SUB1, SUB2；68 行）
%  输出：震源脑区列表、全脑排序表、汇总结果
% =========================================================================

clear; clc; close all;
warning('off', 'all');

%% 1. 参数配置
% ---- Repository paths (no edit needed if you keep this folder layout) ----
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(script_dir);            % repository root
input_file = fullfile(repo_root, 'data', 'Fig4_allROI_GOF_volume.xlsx');
output_dir = fullfile(repo_root, 'output', 'Epicenter_Results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

quantile_thresh = 0.85;      % 分位数阈值（前15%）
topN = 10;                   % 固定前N个脑区
filter_negative = true;      % SUB2只考虑正GOF

output_all_rank = true;

%% 2. 读取数据
if ~exist(input_file, 'file')
    error('找不到文件：%s', input_file);
end
data = readtable(input_file, 'VariableNamingRule', 'preserve');
sub1_gof = data{:, 2};   % SUB1列
sub2_gof = data{:, 3};   % SUB2列
n_roi = length(sub1_gof);
fprintf('成功读取数据，脑区数量：%d\n', n_roi);

%% 3. 筛选函数
function [idx, val, cutoff] = filter_by_quantile(gof, thresh)
    cutoff = quantile(gof, thresh);
    idx = find(gof >= cutoff);
    val = gof(idx);
    [val, order] = sort(val, 'descend');
    idx = idx(order);
end

function [idx, val, cutoff] = filter_by_topN(gof, N, filter_neg)
    if filter_neg
        gof_pos = gof(gof > 0);
        if isempty(gof_pos)
            idx = []; val = []; cutoff = NaN;
            return;
        end
        [sorted_val, sorted_idx_orig] = sort(gof_pos, 'descend');
        orig_idx = find(gof > 0);
        idx = orig_idx(sorted_idx_orig(1:min(N, length(sorted_val))));
        val = sorted_val(1:min(N, length(sorted_val)));
    else
        [sorted_val, sorted_idx] = sort(gof, 'descend');
        idx = sorted_idx(1:min(N, length(sorted_val)));
        val = sorted_val(1:min(N, length(sorted_val)));
    end
    cutoff = val(end);
end

%% 4. 对两个亚型分别处理
subtypes = {'SUB1', 'SUB2'};
gof_cells = {sub1_gof, sub2_gof};
all_results = cell(0, 7);

for s = 1:2
    subtype = subtypes{s};
    gof = gof_cells{s};
    fprintf('\n==================== %s ====================\n', subtype);
    
    % 分位数法
    [idx_q, val_q, cut_q] = filter_by_quantile(gof, quantile_thresh);
    if ~isempty(idx_q)
        fprintf('分位数法 (阈值=%.4f): %d个震源\n', cut_q, length(idx_q));
        fprintf('  脑区: %s\n', mat2str(idx_q'));
        fprintf('  GOF:  %s\n', mat2str(val_q', 4));
        all_results = [all_results; {subtype, 'Quantile', cut_q, mat2str(idx_q'), mat2str(val_q', 4), length(idx_q), ''}];
    else
        all_results = [all_results; {subtype, 'Quantile', cut_q, '', '', 0, '无脑区'}];
    end
    
    % 固定TopN法
    filter_neg_flag = filter_negative && strcmp(subtype, 'SUB2');
    [idx_top, val_top, cut_top] = filter_by_topN(gof, topN, filter_neg_flag);
    if ~isempty(idx_top)
        fprintf('Top%d法 (阈值=%.4f): %d个震源\n', topN, cut_top, length(idx_top));
        fprintf('  脑区: %s\n', mat2str(idx_top'));
        fprintf('  GOF:  %s\n', mat2str(val_top', 4));
        all_results = [all_results; {subtype, 'TopN', cut_top, mat2str(idx_top'), mat2str(val_top', 4), length(idx_top), ''}];
    else
        all_results = [all_results; {subtype, 'TopN', NaN, '', '', 0, '无正GOF'}];
    end
    
    % 输出全脑排序表
    if output_all_rank
        [sorted_val, sorted_idx] = sort(gof, 'descend');
        rank_table = table((1:n_roi)', sorted_idx, sorted_val, ...
            'VariableNames', {'Rank', 'ROI', 'GOF'});
        out_file = fullfile(output_dir, sprintf('%s_AllROI_Ranking.xlsx', subtype));
        writetable(rank_table, out_file);
        fprintf('全脑排序表: %s\n', out_file);
    end
end

%% 5. 保存汇总结果
summary_table = cell2table(all_results, ...
    'VariableNames', {'Subtype', 'Method', 'Cutoff', 'Epicenter_ROIs', 'GOF_Values', 'Count', 'Note'});
out_summary = fullfile(output_dir, 'Epicenter_Summary.xlsx');
writetable(summary_table, out_summary);
fprintf('\n✅ 所有结果已保存至：%s\n', output_dir);