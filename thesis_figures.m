% THESIS_FIGURES  毕业设计论文插图统一重绘脚本
%
% 功能：
%   加载各仿真脚本保存的 .mat 数据，按本科毕设论文排版规范重新绘制
%   所有仿真图，生成可直接插入 Word 的高质量图像文件。
%
% 排版规范：
%   - 中文字体：宋体 (SimSun)，10.5pt (五号)
%   - 英文/数字：Times New Roman，10.5pt
%   - 图例字号：9pt (小五号)
%   - 线宽 1.5pt，标记尺寸 6-8pt
%   - 不同线型/标记确保黑白可辨识
%   - 无标题（标题由 Word 题注系统管理）
%   - 单栏图宽 8.5cm，双栏图宽 17cm
%   - 600dpi PNG + 矢量 PDF 输出
%   - 浅灰网格线 0.5pt
%
% 使用方法：
%   1. 先运行 sim_01 ~ sim_05 各仿真脚本（生成 .mat 数据文件）
%   2. 运行本脚本重新绘制所有论文插图
%
% -----------------------------------------------------------------

clear; clc; close all;

%% =================================================================
%%  全局配置
%% =================================================================
dataFolder   = 'simulation_results';
outputFolder = 'thesis_figures';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% ---- 字体配置 ----
FontCN     = 'SimSun';           % 中文：宋体
FontEN     = 'Times New Roman';  % 英文/数字
FontSize   = 10.5;               % 五号字
LegendSize = 9;                  % 小五号

% ---- 线条配置 ----
LW  = 1.5;   % 线宽 (pt)
MS  = 7;     % 标记尺寸 (pt)
GridLW = 0.5; % 网格线宽

% ---- 图像尺寸 (cm) ----
SingleW = 8.5;   % 单栏
DoubleW = 17;    % 双栏
FigH    = 6.5;   % 单图高度

% ---- 配色方案（学术规范，黑白可辨识）----
blue   = [0.0  0.45 0.75];
red    = [0.75 0.0  0.0 ];
green  = [0.0  0.6  0.3 ];
black  = [0.05 0.05 0.05];
gray   = [0.5  0.5  0.5 ];

%% =================================================================
%%  辅助函数：统一坐标轴格式
%% =================================================================
setupAxes = @(ax) set(ax, ...
    'FontName', FontEN, 'FontSize', FontSize, ...
    'TickDir', 'in', 'TickLength', [0.02 0.02], ...
    'Box', 'on', 'LineWidth', 0.5, ...
    'GridLineStyle', '-', 'GridAlpha', 0.15, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridLineStyle', ':', 'MinorGridAlpha', 0.08, ...
    'XMinorTick', 'on', 'YMinorTick', 'on');

%% =================================================================
%%  辅助函数：设置标签字体（中英文混排）
%% =================================================================
setLabel = @(h, str) set(h, 'String', str, ...
    'FontName', FontCN, 'FontSize', FontSize, 'FontWeight', 'normal');

%% =================================================================
%%  辅助函数：导出图像
%% =================================================================
function exportFig(fig, folder, name, widthCM, heightCM)
    set(fig, 'PaperUnits', 'centimeters', ...
             'PaperSize', [widthCM heightCM], ...
             'PaperPosition', [0 0 widthCM heightCM]);
    exportgraphics(fig, fullfile(folder, [name, '.png']), ...
        'Resolution', 600, 'BackgroundColor', 'white');
    exportgraphics(fig, fullfile(folder, [name, '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('  => %s.png / .pdf\n', name);
end


%% #####################################################################
%%  图1：V2V SINR 累积分布函数 (sim_01)
%% #####################################################################
fprintf('绘制图1：V2V SINR CDF ...\n');
d1 = load(fullfile(dataFolder, 'sim_01_data.mat'));

% 排序计算经验 CDF
[sR,  ~] = sort(10*log10(d1.sinr_robust));
cdfR = (1:length(sR))' / length(sR);
[sNR, ~] = sort(10*log10(d1.sinr_nonrobust));
cdfNR = (1:length(sNR))' / length(sNR);

% 均匀抽样控制阶梯密度
stepR  = unique(round(linspace(1, length(sR),  3000)));
stepNR = unique(round(linspace(1, length(sNR), 3000)));

fig1 = figure('Visible', 'off', 'Color', 'w');
ax1 = axes(fig1);
setupAxes(ax1);
grid(ax1, 'on'); hold(ax1, 'on');

stairs(ax1, sR(stepR),   cdfR(stepR),   '-',  'Color', blue, 'LineWidth', LW);
stairs(ax1, sNR(stepNR), cdfNR(stepNR), '--', 'Color', red,  'LineWidth', LW);

% SINR 阈值参考线
plot(ax1, [5 5], [0 1], ':', 'Color', gray, 'LineWidth', 1);

% 中断概率标记
plot(ax1, 5, d1.outage_robust,    'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', blue, 'MarkerEdgeColor', 'none');
plot(ax1, 5, d1.outage_nonrobust, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', red,  'MarkerEdgeColor', 'none');

% 标注文字
text(ax1, 6, 0.15, sprintf('%.2f%%', d1.outage_robust*100), ...
    'FontName', FontEN, 'FontSize', LegendSize, 'Color', blue);
text(ax1, 6, 0.80, sprintf('%.1f%%', d1.outage_nonrobust*100), ...
    'FontName', FontEN, 'FontSize', LegendSize, 'Color', red);

setLabel(xlabel(ax1, ''), 'V2V SINR (dB)');
setLabel(ylabel(ax1, ''), {'\fontname{SimSun}累积分布函数', ' (CDF)'});
ylim(ax1, [0 1.05]);
xmin = min(min(sR), min(sNR));
xlim(ax1, [floor(xmin-2), max(max(sR), max(sNR))+3]);

lg1 = legend(ax1, {'\fontname{SimSun}鲁棒算法', '\fontname{SimSun}非鲁棒算法'}, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'southeast', 'Box', 'off');
hold(ax1, 'off');

exportFig(fig1, outputFolder, 'fig_01_V2V_SINR_CDF', SingleW, FigH);


%% #####################################################################
%%  图2：有效V2I吞吐量 vs CSI反馈周期 (sim_02)
%% #####################################################################
fprintf('绘制图2：V2I吞吐量 vs CSI周期 ...\n');
d2 = load(fullfile(dataFolder, 'sim_02_data.mat'));

colors3 = {blue, black, green};
markers3 = {'o', 's', '^'};
lineR  = {'-',  '-',  '-'};
lineNR = {'--', '--', '--'};

fig2 = figure('Visible', 'off', 'Color', 'w');
ax2 = axes(fig2);
setupAxes(ax2);
grid(ax2, 'on'); hold(ax2, 'on');

h_leg = []; leg_str = {};

% 鲁棒算法（实线 + 填充标记）
for vi = 1 : length(d2.v_list)
    c = colors3{vi}; mk = markers3{vi};
    h = plot(ax2, d2.T_list, d2.sumRate_robust(vi,:), ...
        [mk, lineR{vi}], 'LineWidth', LW, 'MarkerSize', MS, ...
        'MarkerFaceColor', c, 'Color', c);
    h_leg(end+1) = h;
    leg_str{end+1} = sprintf('\\fontname{SimSun}鲁棒, v=%d km/h', d2.v_list(vi));
end

% 非鲁棒算法（虚线 + 空心标记）
for vi = 1 : length(d2.v_list)
    c = colors3{vi}; mk = markers3{vi};
    h = plot(ax2, d2.T_list, d2.sumRate_nonrobust(vi,:), ...
        [mk, lineNR{vi}], 'LineWidth', LW, 'MarkerSize', MS, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', c, 'Color', c);
    h_leg(end+1) = h;
    leg_str{end+1} = sprintf('\\fontname{SimSun}非鲁棒, v=%d km/h', d2.v_list(vi));
end

setLabel(xlabel(ax2, ''), 'CSI\fontname{SimSun}反馈周期\fontname{Times New Roman} T (ms)');
setLabel(ylabel(ax2, ''), '\fontname{SimSun}有效\fontname{Times New Roman}V2I\fontname{SimSun}总吞吐量\fontname{Times New Roman} (bps/Hz)');
xlim(ax2, [0 5]);
ylim(ax2, [0 max(d2.sumRate_robust(:))*1.15]);

legend(ax2, h_leg, leg_str, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'northeast', 'Box', 'off');
hold(ax2, 'off');

exportFig(fig2, outputFolder, 'fig_02_V2I_Rate_vs_T', SingleW, FigH);


%% #####################################################################
%%  图3：V2V中断概率和V2I容量 vs 密度 (sim_03)  — 双栏子图
%% #####################################################################
fprintf('绘制图3：密度对V2V中断和V2I容量的影响 ...\n');
d3 = load(fullfile(dataFolder, 'sim_03_data.mat'));

fig3 = figure('Visible', 'off', 'Color', 'w');

% ---- 子图(a)：V2V 中断概率 ----
ax3a = subplot(1, 2, 1, 'Parent', fig3);
setupAxes(ax3a);
grid(ax3a, 'on'); hold(ax3a, 'on');

h_r = semilogy(ax3a, d3.N_list, d3.P_outage_robust, 'o-', ...
    'LineWidth', LW, 'MarkerSize', MS+1, ...
    'MarkerFaceColor', blue, 'Color', blue);
h_nr = semilogy(ax3a, d3.N_list, d3.P_outage_nonrobust, 's--', ...
    'LineWidth', LW, 'MarkerSize', MS+1, ...
    'MarkerFaceColor', red, 'Color', red);

setLabel(xlabel(ax3a, ''), '\fontname{SimSun}车辆密度\fontname{Times New Roman} N (\fontname{SimSun}辆\fontname{Times New Roman})');
setLabel(ylabel(ax3a, ''), 'V2V\fontname{SimSun}链路中断概率');
set(ax3a, 'YScale', 'log');
ylim(ax3a, [1e-3 1]);
xlim(ax3a, [5 50]);

text(ax3a, 0.05, 0.05, '(a)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

legend(ax3a, [h_r, h_nr], ...
    {'\fontname{SimSun}鲁棒算法', '\fontname{SimSun}非鲁棒算法'}, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'southeast', 'Box', 'off');

% ---- 子图(b)：总有效V2I容量 ----
ax3b = subplot(1, 2, 2, 'Parent', fig3);
setupAxes(ax3b);
grid(ax3b, 'on'); hold(ax3b, 'on');

plot(ax3b, d3.N_list, d3.TotalCap_robust, '^-', ...
    'LineWidth', LW, 'MarkerSize', MS+1, ...
    'MarkerFaceColor', blue, 'Color', blue);
plot(ax3b, d3.N_list, d3.TotalCap_nonrobust, 'd--', ...
    'LineWidth', LW, 'MarkerSize', MS+1, ...
    'MarkerFaceColor', red, 'Color', red);

setLabel(xlabel(ax3b, ''), '\fontname{SimSun}车辆密度\fontname{Times New Roman} N (\fontname{SimSun}辆\fontname{Times New Roman})');
setLabel(ylabel(ax3b, ''), '\fontname{SimSun}总有效\fontname{Times New Roman}V2I\fontname{SimSun}容量\fontname{Times New Roman} (bps/Hz)');
ylim(ax3b, [0 max(max(d3.TotalCap_robust), max(d3.TotalCap_nonrobust))*1.2]);
xlim(ax3b, [5 50]);

text(ax3b, 0.05, 0.05, '(b)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

legend(ax3b, {'\fontname{SimSun}鲁棒算法', '\fontname{SimSun}非鲁棒算法'}, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'best', 'Box', 'off');
hold(ax3a, 'off'); hold(ax3b, 'off');

exportFig(fig3, outputFolder, 'fig_03_Density_Outage_Cap', DoubleW, FigH);


%% #####################################################################
%%  图4：收敛性与计算复杂度 (sim_04)  — 双栏子图
%% #####################################################################
fprintf('绘制图4：算法收敛性与复杂度 ...\n');
d4 = load(fullfile(dataFolder, 'sim_04_data.mat'));

fig4 = figure('Visible', 'off', 'Color', 'w');

% ---- 子图(a)：收敛曲线 ----
ax4a = subplot(1, 2, 1, 'Parent', fig4);
setupAxes(ax4a);
grid(ax4a, 'on'); hold(ax4a, 'on');

iter_axis = 1 : d4.max_iters;
plot(ax4a, iter_axis, d4.mean_obj_hist, '-', ...
    'Color', blue, 'LineWidth', LW);
plot(ax4a, [1 d4.max_iters], [d4.mean_obj_nr d4.mean_obj_nr], '--', ...
    'Color', red, 'LineWidth', LW);
plot(ax4a, [1 d4.max_iters], [d4.mean_obj_final_r d4.mean_obj_final_r], ':', ...
    'Color', green, 'LineWidth', LW);
plot(ax4a, d4.avg_iters, d4.mean_obj_final_r, '^', ...
    'MarkerSize', MS+1, 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', blue, 'LineWidth', 1.2);

setLabel(xlabel(ax4a, ''), '\fontname{SimSun}迭代次数');
setLabel(ylabel(ax4a, ''), 'V2I\fontname{SimSun}容量\fontname{Times New Roman} (bps/Hz)');
xlim(ax4a, [0 d4.max_iters+2]);
ylim(ax4a, [0 max(d4.mean_obj_nr, d4.mean_obj_final_r)*1.25]);

text(ax4a, 0.05, 0.05, '(a)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

legend(ax4a, ...
    {'\fontname{SimSun}鲁棒算法（二分搜索）', ...
     '\fontname{SimSun}非鲁棒算法（闭式解）', ...
     sprintf('\\fontname{SimSun}鲁棒收敛值 (iter=%d)', d4.avg_iters)}, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'southeast', 'Box', 'off');

% ---- 子图(b)：平均迭代次数 vs 密度 ----
ax4b = subplot(1, 2, 2, 'Parent', fig4);
setupAxes(ax4b);
grid(ax4b, 'on'); hold(ax4b, 'on');

plot(ax4b, d4.N_list, d4.avg_iter_robust, '^-', ...
    'LineWidth', LW, 'MarkerSize', MS+2, ...
    'MarkerFaceColor', blue, 'Color', blue);
plot(ax4b, d4.N_list, d4.avg_iter_nonrobust, 'd--', ...
    'LineWidth', LW, 'MarkerSize', MS+2, ...
    'MarkerFaceColor', red, 'Color', red);

setLabel(xlabel(ax4b, ''), '\fontname{SimSun}车辆密度\fontname{Times New Roman} N');
setLabel(ylabel(ax4b, ''), '\fontname{SimSun}平均迭代次数');
xlim(ax4b, [5 50]);
ylim(ax4b, [-2 max(d4.avg_iter_robust)*1.3]);

text(ax4b, 0.05, 0.05, '(b)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

legend(ax4b, {'\fontname{SimSun}鲁棒算法', '\fontname{SimSun}非鲁棒算法'}, ...
    'FontName', FontCN, 'FontSize', LegendSize, ...
    'Location', 'northwest', 'Box', 'off');
hold(ax4a, 'off'); hold(ax4b, 'off');

exportFig(fig4, outputFolder, 'fig_04_Convergence_Complexity', DoubleW, FigH);


%% #####################################################################
%%  图5：SINR阈值敏感性分析 (sim_05)  — 双栏子图
%% #####################################################################
fprintf('绘制图5：SINR阈值敏感性 ...\n');
d5 = load(fullfile(dataFolder, 'sim_05_data.mat'));

% 配色：蓝色系（鲁棒，从浅到深）、红色系（非鲁棒，从浅到深）
blues = [0.7 0.85 1.0; 0.5 0.7 0.9; 0.2 0.5 0.8; 0.1 0.35 0.65; 0.0 0.2 0.5];
reds  = [1.0 0.7 0.7; 0.9 0.5 0.5; 0.75 0.25 0.25; 0.6 0.1 0.1; 0.45 0.0 0.0];
mk5 = {'o', 's', '^', 'd', 'v'};
ls5_r  = {'-',  '-',  '-',  '-',  '-'};
ls5_nr = {'--', '--', '--', '--', '--'};

fig5 = figure('Visible', 'off', 'Color', 'w');

% ---- 子图(a)：鲁棒算法 ----
ax5a = subplot(1, 2, 1, 'Parent', fig5);
setupAxes(ax5a);
grid(ax5a, 'on'); hold(ax5a, 'on');

for ti = 1 : length(d5.T_list)
    c = blues(ti,:);
    semilogy(ax5a, d5.gamma_th_dB_list, d5.outage_prob_robust_raw(:,ti), ...
        [mk5{ti}, ls5_r{ti}], 'LineWidth', LW, 'MarkerSize', MS, ...
        'MarkerFaceColor', c, 'Color', c);
end

% 垂直参考线
plot(ax5a, [5 5], [1e-3 1], ':', 'Color', gray, 'LineWidth', 1);

setLabel(xlabel(ax5a, ''), 'SINR\fontname{SimSun}阈值\fontname{Times New Roman} \gamma_{th} (dB)');
setLabel(ylabel(ax5a, ''), 'V2V\fontname{SimSun}中断概率\fontname{Times New Roman} P_{out}');
set(ax5a, 'YScale', 'log');
ylim(ax5a, [1e-3 1]); xlim(ax5a, [-6 16]);

text(ax5a, 0.05, 0.05, '(a)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

leg5a = arrayfun(@(t) sprintf('T=%.1f ms', t), d5.T_list, 'UniformOutput', false);
legend(ax5a, leg5a, 'FontName', FontEN, 'FontSize', LegendSize, ...
    'Location', 'southwest', 'Box', 'off');

% ---- 子图(b)：非鲁棒算法 ----
ax5b = subplot(1, 2, 2, 'Parent', fig5);
setupAxes(ax5b);
grid(ax5b, 'on'); hold(ax5b, 'on');

for ti = 1 : length(d5.T_list)
    c = reds(ti,:);
    semilogy(ax5b, d5.gamma_th_dB_list, d5.outage_prob_nonrobust_raw(:,ti), ...
        [mk5{ti}, ls5_nr{ti}], 'LineWidth', LW, 'MarkerSize', MS, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', c, 'Color', c);
end

plot(ax5b, [5 5], [1e-3 1], ':', 'Color', gray, 'LineWidth', 1);

setLabel(xlabel(ax5b, ''), 'SINR\fontname{SimSun}阈值\fontname{Times New Roman} \gamma_{th} (dB)');
setLabel(ylabel(ax5b, ''), 'V2V\fontname{SimSun}中断概率\fontname{Times New Roman} P_{out}');
set(ax5b, 'YScale', 'log');
ylim(ax5b, [1e-3 1]); xlim(ax5b, [-6 16]);

text(ax5b, 0.05, 0.05, '(b)', 'Units', 'normalized', ...
    'FontName', FontEN, 'FontSize', FontSize, 'FontWeight', 'bold');

leg5b = arrayfun(@(t) sprintf('T=%.1f ms', t), d5.T_list, 'UniformOutput', false);
legend(ax5b, leg5b, 'FontName', FontEN, 'FontSize', LegendSize, ...
    'Location', 'northwest', 'Box', 'off');
hold(ax5a, 'off'); hold(ax5b, 'off');

exportFig(fig5, outputFolder, 'fig_05_SINR_Threshold', DoubleW, FigH);


%% =================================================================
%%  完成
%% =================================================================
fprintf('\n===== 所有论文插图已生成 =====\n');
fprintf('输出目录: %s/\n', outputFolder);
fprintf('共 5 张图，每张 PNG (600dpi) + PDF (矢量)\n');
close all;
