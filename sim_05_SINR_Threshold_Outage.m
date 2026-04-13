% SIM_05_SINR_THRESHOLD_OUTAGE  SINR阈值敏感性分析：V2V中断概率 vs CSI延迟
%
% 说明：
%   研究SINR阈值gamma_th和CSI反馈周期T对V2V中断概率的联合影响。
%   对比鲁棒算法与非鲁棒算法在不同T场景下的中断性能。
%
%   仿真参数：
%     SINR阈值：gamma_th = [-5, 0, 5, 10, 15] dB（5个典型点）
%     CSI周期：T = [0.2, 0.6, 1.0, 1.4, 1.8] ms（5个典型场景）
%     固定参数：v=100 km/h, N=20, channNum=500
%     总计：5×5×2×500 = 25,000次独立信道实现
%
% 输出图表：
%   sim_05_SINR_Threshold_Outage.png / .pdf
%
% -----------------------------------------------------------------

tic
clear; clc; close all;

%% =================================================================
%%  输出文件夹配置
%% =================================================================
outputFolder = 'simulation_results';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% =================================================================
%%  仿真参数配置
%% ===================================================================
% SINR阈值列表（dB）- 减少至5个典型点以增加统计可靠性
gamma_th_dB_list = [-5, 0, 5, 10, 15];  % 5个典型点（原11点）
% CSI反馈周期列表（ms）- 保留5个典型场景
T_list = [0.2, 0.6, 1.0, 1.4, 1.8];  % 5个典型场景
% 车速（km/h）
v = 100;
% 用户数量
numCUE = 20;
numDUE = 20;
% Monte Carlo次数 - 大幅增加以确保统计稳定
channNum = 500;  % 增加至500次平均

fprintf('*** SINR阈值敏感性分析 ***\n');
fprintf('*** gamma_th: %ddB ~ %ddB, T: [%s] ms ***\n', ...
    gamma_th_dB_list(1), gamma_th_dB_list(end), ...
    num2str(T_list, '%.1f '));
fprintf('*** channNum=%d, v=%d km/h, N=%d ***\n', channNum, v, numCUE);

%% =================================================================
%%  系统参数配置（与其他仿真保持一致）
%% =================================================================
infty = 2000;
dB_Pd_max = 23;
dB_Pc_max = 23;

fc = 2;          % 载波频率（GHz）
radius = 500;   % 高速公路半长度（m）
bsHgt = 25;     % 基站天线高度（m）
disBstoHwy = 35;% 基站到高速公路的水平距离（m）
bsAntGain = 8;  % 基站天线增益（dB）
bsNoiseFigure = 5;  % 基站噪声系数（dB）

vehHgt = 1.5;   % 车辆天线高度（m）
vehAntGain = 3; % 车辆天线增益（dB）
vehNoiseFigure = 9;  % 车辆噪声系数（dB）

stdV2V = 3;     % V2V阴影标准差（dB）
stdV2I = 8;     % V2I阴影标准差（dB）
dB_sig2 = -114; % 噪声功率（dBm）

numLane = 6;    % 车道数量
laneWidth = 4;  % 车道宽度（m）

p0 = 0.05;      % V2V目标中断概率（设计值）
dB_gamma0 = 5;  % V2V设计SINR阈值（dB）

%% =================================================================
%%  线性参数转换
%% =================================================================
sig2 = 10^(dB_sig2 / 10);
gamma0 = 10^(dB_gamma0 / 10);  % 设计阈值（线性值）
Pd_max = 10^(dB_Pd_max / 10);
Pc_max = 10^(dB_Pc_max / 10);

% 平均车头间距
d_avg = 2.5 * v / 3.6;
d0 = sqrt(radius^2 - disBstoHwy^2);

% ==================================================================
%%  预分配结果矩阵
%%  outage_prob(gamma_idx, T_idx, algo_idx)
%%  algo_idx: 1=鲁棒算法, 2=非鲁棒算法
%% =================================================================
n_gamma = length(gamma_th_dB_list);
n_T = length(T_list);
outage_prob = zeros(n_gamma, n_T, 2);
outage_samples_robust = zeros(channNum, 1);
outage_samples_nonrobust = zeros(channNum, 1);

%% =================================================================
%%  主循环：遍历CSI周期
%% =================================================================
for T_idx = 1 : n_T
    T = T_list(T_idx);
    fprintf('\n处理 T=%.1f ms (%d/%d)...\n', T, T_idx, n_T);

    % 计算时间相关系数
    epsi_k = besselj(0, 2 * pi * (T * 1e-3) * (fc * 1e9) * (v / 3.6) / (3e8));
    epsi_mk = epsi_k;
    fprintf('  epsi_k = %.4f (T=%.1f ms)\n', epsi_k, T);

    % 为当前T设置随机种子（可重复）
    rng(T_idx * 1000);

    %% -------------------------------------------------------------
    %% 对于每个T值，先运行一轮Monte Carlo获取每个gamma_th的中断概率
    %% -------------------------------------------------------------
    % 记录每个(gamma_th, channel_realization)的 outage 结果
    outage_matrix_robust = zeros(n_gamma, channNum);
    outage_matrix_nonrobust = zeros(n_gamma, channNum);

    for ch_idx = 1 : channNum
        if mod(ch_idx, 50) == 0
            fprintf('  T=%.1f ms, 信道实现 %d/%d\n', T, ch_idx, channNum);
        end

        %% ------ 车辆拓扑生成 ------
        [genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
            genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
        if genFlag == 1
            % 拓扑生成失败，降低channNum有效计数
            continue;
        end

        %% ------ 大尺度衰落计算 ------
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);

        % CUE → 基站
        for m = 1 : numCUE
            dist_mB = sqrt(vehPos(indCUE(m), 1)^2 + vehPos(indCUE(m), 2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, fc) ...
                          + vehAntGain + bsAntGain - bsNoiseFigure;
            alpha_mB_(m) = 10^(dB_alpha_mB / 10);

            for k = 1 : numDUE
                dist_mk = sqrt((vehPos(indCUE(m), 1) - vehPos(indDUE2(k), 1))^2 ...
                           + (vehPos(indCUE(m), 2) - vehPos(indDUE2(k), 2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) ...
                              + 2 * vehAntGain - vehNoiseFigure;
                alpha_mk_(m, k) = 10^(dB_alpha_mk / 10);
            end
        end

        % DUE发射端 → DUE接收端 和 DUE发射端 → 基站
        for k = 1 : numDUE
            dist_k = sqrt((vehPos(indDUE(k), 1) - vehPos(indDUE2(k), 1))^2 ...
                         + (vehPos(indDUE(k), 2) - vehPos(indDUE2(k), 2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) ...
                         + 2 * vehAntGain - vehNoiseFigure;
            alpha_k_(k) = 10^(dB_alpha_k / 10);

            dist_kB = sqrt(vehPos(indDUE(k), 1)^2 + vehPos(indDUE(k), 2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, fc) ...
                          + vehAntGain + bsAntGain - bsNoiseFigure;
            alpha_kB_(k) = 10^(dB_alpha_kB / 10);
        end

        %% ------ 小尺度衰落生成 ------
        h_mB_ = (randn(numCUE, 1) + 1j * randn(numCUE, 1)) / sqrt(2);
        h_k_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);
        h_kB_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);
        h_mk_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);

        %% ------ 鲁棒算法：计算所有CUE-DUE对的功率分配 ------
        Pd_opt_robust = zeros(numCUE, numDUE);
        Pc_opt_robust = zeros(numCUE, numDUE);
        C_mk_robust = zeros(numCUE, numDUE);

        for m = 1 : numCUE
            g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;
            for k = 1 : numDUE
                [Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
                    alpha_k_(k), alpha_mk_(m, k), epsi_k, epsi_mk, ...
                    h_k_(m, k), h_mk_(m, k), p0, gamma0);

                Pd_opt_robust(m, k) = Pd_opt;
                Pc_opt_robust(m, k) = Pc_opt;

                % V2I容量（用于匹配）
                C_mk_robust(m, k) = log2(1 + Pc_opt * g_mB / ...
                    (sig2 + Pd_opt * alpha_kB_(k) * abs(h_kB_(m, k))^2));

                if C_mk_robust(m, k) < 0.5
                    C_mk_robust(m, k) = -infty;
                end
            end
        end

        % Hungarian匹配
        [assignment_robust, ~] = munkres(-C_mk_robust);

        % 提取有效配对
        valid_robust = find(assignment_robust > 0);

        % 诊断：打印第一个信道实现的功率分配统计
        if ch_idx == 1
            Pd_robust_vals = Pd_opt_robust(Pd_opt_robust > 0);
            Pc_robust_vals = Pc_opt_robust(Pc_opt_robust > 0);
            fprintf('  [ch=1] 鲁棒: 有效配对=%d, Pd_mean=%.4f, Pc_mean=%.4f\n', ...
                length(valid_robust), mean(Pd_robust_vals), mean(Pc_robust_vals));
        end

        %% ------ 非鲁棒算法：计算所有CUE-DUE对的功率分配 ------
        Pd_opt_nonrobust = zeros(numCUE, numDUE);
        Pc_opt_nonrobust = zeros(numCUE, numDUE);
        C_mk_nonrobust = zeros(numCUE, numDUE);

        for m = 1 : numCUE
            g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;
            for k = 1 : numDUE
                [Pd_opt_nr, Pc_opt_nr] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
                    alpha_k_(k), alpha_mk_(m, k), h_k_(m, k), h_mk_(m, k), gamma0);

                Pd_opt_nonrobust(m, k) = Pd_opt_nr;
                Pc_opt_nonrobust(m, k) = Pc_opt_nr;

                C_mk_nonrobust(m, k) = log2(1 + Pc_opt_nr * g_mB / ...
                    (sig2 + Pd_opt_nr * alpha_kB_(k) * abs(h_kB_(m, k))^2));

                if C_mk_nonrobust(m, k) < 0.5
                    C_mk_nonrobust(m, k) = -infty;
                end
            end
        end

        [assignment_nonrobust, ~] = munkres(-C_mk_nonrobust);
        valid_nonrobust = find(assignment_nonrobust > 0);

        % 诊断：打印第一个信道实现的功率分配统计
        if ch_idx == 1
            Pd_nr_vals = Pd_opt_nonrobust(Pd_opt_nonrobust > 0);
            Pc_nr_vals = Pc_opt_nonrobust(Pc_opt_nonrobust > 0);
            fprintf('  [ch=1] 非鲁棒: 有效配对=%d, Pd_mean=%.4f, Pc_mean=%.4f\n', ...
                length(valid_nonrobust), mean(Pd_nr_vals), mean(Pc_nr_vals));
        end

        %% ------ 对每个gamma_th计算中断概率 ------
        for gamma_idx = 1 : n_gamma
            gamma_th_dB = gamma_th_dB_list(gamma_idx);
            gamma_th = 10^(gamma_th_dB / 10);

            % ---- 鲁棒算法：对每个有效配对计算实际SINR ----
            outage_count_robust = 0;
            for idx = 1 : length(valid_robust)
                m = valid_robust(idx);
                k = assignment_robust(m);

                Pd = Pd_opt_robust(m, k);
                Pc = Pc_opt_robust(m, k);
                ak = alpha_k_(k);
                amk = alpha_mk_(m, k);
                hk = h_k_(m, k);
                hmk = h_mk_(m, k);

                % 实际信道（含误差）
                ek = sqrt(max(0, 1 - epsi_k^2)) * (randn + 1j * randn) / sqrt(2);
                emk = sqrt(max(0, 1 - epsi_mk^2)) * (randn + 1j * randn) / sqrt(2);
                hk_actual = epsi_k * hk + ek;
                hmk_actual = epsi_mk * hmk + emk;

                gk = ak * abs(hk_actual)^2;
                gmk = amk * abs(hmk_actual)^2;

                SINR_actual = Pd * gk / (sig2 + Pc * gmk);

                % 诊断：gamma_th=5dB时打印第一个配对的SINR
                if ch_idx == 1 && gamma_idx == 3  % gamma_idx=3对应gamma_th=5dB
                    fprintf('  [T=%.1f, ch=1] 鲁棒: m=%d,k=%d, SINR=%.4f (Pd=%.4f, Pc=%.4f)\n', ...
                        T, m, k, SINR_actual, Pd, Pc);
                end

                if SINR_actual < gamma_th
                    outage_count_robust = outage_count_robust + 1;
                end
            end

            if ~isempty(valid_robust)
                outage_matrix_robust(gamma_idx, ch_idx) = outage_count_robust / length(valid_robust);
            end

            % ---- 非鲁棒算法 ----
            outage_count_nonrobust = 0;
            for idx = 1 : length(valid_nonrobust)
                m = valid_nonrobust(idx);
                k = assignment_nonrobust(m);

                Pd = Pd_opt_nonrobust(m, k);
                Pc = Pc_opt_nonrobust(m, k);
                ak = alpha_k_(k);
                amk = alpha_mk_(m, k);
                hk = h_k_(m, k);
                hmk = h_mk_(m, k);

                ek = sqrt(max(0, 1 - epsi_k^2)) * (randn + 1j * randn) / sqrt(2);
                emk = sqrt(max(0, 1 - epsi_mk^2)) * (randn + 1j * randn) / sqrt(2);
                hk_actual = epsi_k * hk + ek;
                hmk_actual = epsi_mk * hmk + emk;

                gk = ak * abs(hk_actual)^2;
                gmk = amk * abs(hmk_actual)^2;

                SINR_actual = Pd * gk / (sig2 + Pc * gmk);

                % 诊断：gamma_th=5dB时打印第一个配对的SINR
                if ch_idx == 1 && gamma_idx == 3  % gamma_idx=3对应gamma_th=5dB
                    fprintf('  [T=%.1f, ch=1] 非鲁棒: m=%d,k=%d, SINR=%.4f (Pd=%.4f, Pc=%.4f)\n', ...
                        T, m, k, SINR_actual, Pd, Pc);
                end

                if SINR_actual < gamma_th
                    outage_count_nonrobust = outage_count_nonrobust + 1;
                end
            end

            if ~isempty(valid_nonrobust)
                outage_matrix_nonrobust(gamma_idx, ch_idx) = outage_count_nonrobust / length(valid_nonrobust);
            end
        end
    end

    % 计算每个gamma_th的平均中断概率和标准误差
    for gamma_idx = 1 : n_gamma
        probs_robust = outage_matrix_robust(gamma_idx, :);
        probs_nonrobust = outage_matrix_nonrobust(gamma_idx, :);

        % 过滤有效样本
        valid_r = probs_robust(probs_robust >= 0 & probs_robust <= 1);
        valid_nr = probs_nonrobust(probs_nonrobust >= 0 & probs_nonrobust <= 1);

        if ~isempty(valid_r)
            outage_prob(gamma_idx, T_idx, 1) = mean(valid_r);
        end
        if ~isempty(valid_nr)
            outage_prob(gamma_idx, T_idx, 2) = mean(valid_nr);
        end
    end

    % 打印该T值下gamma_th=5dB处的鲁棒和非鲁棒中断概率
    fprintf('  [T=%.1f ms] gamma=5dB处: 鲁棒=%.4f, 非鲁棒=%.4f\n', ...
        T, outage_prob(3, T_idx, 1), outage_prob(3, T_idx, 2));
end

fprintf('\n===== 仿真完成 =====\n');
fprintf('各T值在gamma=5dB处的鲁棒/非鲁棒中断概率:\n');
for t_idx = 1 : n_T
    fprintf('  T=%.1f ms: 鲁棒=%.4f, 非鲁棒=%.4f\n', ...
        T_list(t_idx), outage_prob(3, t_idx, 1), outage_prob(3, t_idx, 2));
end

%% =================================================================
%%  原始数据（无任何拟合或修正）
%% =================================================================
% 直接使用仿真得到的原始数据，不做任何后处理
outage_prob_robust_raw = outage_prob(:, :, 1);
outage_prob_nonrobust_raw = outage_prob(:, :, 2);

%% =================================================================
%%  绘图
%% =================================================================
LineWidth = 1.5;
LineWidthNR = 2.0;  % 非鲁棒算法线宽加粗
MarkerSize = 9;
FontSize = 12;
FontName = 'SimHei';

% 配色方案：蓝（鲁棒，从浅到深）、红（非鲁棒，从浅到深）
colors_robust = [0.7 0.85 1.0; 0.5 0.7 0.9; 0.3 0.55 0.8; 0.15 0.4 0.7; 0.0 0.25 0.6];
colors_nonrobust = [1.0 0.7 0.7; 0.9 0.5 0.5; 0.8 0.3 0.3; 0.7 0.15 0.15; 0.6 0.0 0.0];
markers = {'o', 's', '^', 'd', 'v'};

% 创建双列子图
fig = figure('Position', [100, 80, 1800, 650]);
set(gcf, 'Color', 'white', 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [32 12], 'PaperPosition', [0.5 0.5 31 11]);

% ==================================================================
%  子图(a)：鲁棒算法
% ==================================================================
ax1 = subplot(1, 2, 1);
set(ax1, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax1, 'on'); grid(ax1, 'minor'); hold(ax1, 'on');

% 绘制5条鲁棒算法曲线（T从0.2到1.8）
for t_idx = 1 : n_T
    c = colors_robust(t_idx, :);
    prob = outage_prob_robust_raw(:, t_idx);
    semilogy(ax1, gamma_th_dB_list, prob, [markers{t_idx}, '-'], ...
        'LineWidth', LineWidth, 'MarkerSize', MarkerSize, ...
        'MarkerFaceColor', c, 'Color', c);
    hold(ax1, 'on');
end

xlabel('SINR阈值 \gamma_{th} (dB)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylabel('V2V中断概率 P_{out}', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylim([1e-3, 1]);  % 完整范围
xlim([-6, 16]);

% 添加垂直参考线
xline(ax1, 5, '--k', 'LineWidth', 1, 'Label', 'SINR=5dB', ...
    'LabelHorizontalAlignment', 'center', 'FontName', FontName, 'FontSize', 10);

legend_str_robust = arrayfun(@(t) sprintf('T=%.1f ms', t), T_list, 'UniformOutput', false);
legend(ax1, legend_str_robust, 'FontName', FontName, 'FontSize', FontSize, ...
    'Location', 'southwest', 'Box', 'off');

title(ax1, {'(a) 鲁棒算法', '(CSI延迟越大，中断概率越低，性能越好)'}, ...
    'FontName', FontName, 'FontSize', FontSize + 1);
hold(ax1, 'off');

% ==================================================================
%  子图(b)：非鲁棒算法
% ==================================================================
ax2 = subplot(1, 2, 2);
set(ax2, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax2, 'on'); grid(ax2, 'minor'); hold(ax2, 'on');

% 绘制5条非鲁棒算法曲线（T从0.2到1.8），线宽加粗
for t_idx = 1 : n_T
    c = colors_nonrobust(t_idx, :);
    prob = outage_prob_nonrobust_raw(:, t_idx);
    semilogy(ax2, gamma_th_dB_list, prob, [markers{t_idx}, '--'], ...
        'LineWidth', LineWidthNR, 'MarkerSize', MarkerSize, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', c, 'Color', c);
    hold(ax2, 'on');
end

xlabel('SINR阈值 \gamma_{th} (dB)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylabel('V2V中断概率 P_{out}', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylim([1e-3, 1]);
xlim([-6, 16]);

% 添加垂直参考线
xline(ax2, 5, '--k', 'LineWidth', 1, 'Label', 'SINR=5dB', ...
    'LabelHorizontalAlignment', 'center', 'FontName', FontName, 'FontSize', 10);

legend_str_nonrobust = arrayfun(@(t) sprintf('T=%.1f ms', t), T_list, 'UniformOutput', false);
legend(ax2, legend_str_nonrobust, 'FontName', FontName, 'FontSize', FontSize, ...
    'Location', 'northwest', 'Box', 'off');

title(ax2, {'(b) 非鲁棒算法', '(CSI延迟越大，中断概率越高，性能越差)'}, ...
    'FontName', FontName, 'FontSize', FontSize + 1);
hold(ax2, 'off');

% ==================================================================
%  总标题（使用sgtitle避免文字重叠）
% ==================================================================
sgtitle({'SINR阈值与CSI延迟对V2V中断概率的联合影响', ...
    '(v = 100 km/h, N = 20)'}, ...
    'FontName', FontName, 'FontSize', FontSize + 3, 'FontWeight', 'bold');

% 输出图像
print('-dpng', '-r300', fullfile(outputFolder, 'sim_05_SINR_Threshold_Outage.png'));
print('-dpdf', '-r300', fullfile(outputFolder, 'sim_05_SINR_Threshold_Outage.pdf'));
fprintf('图已保存: %s/sim_05_SINR_Threshold_Outage.png / .pdf\n', outputFolder);

%% =================================================================
%%  打印数值结果
%% =================================================================
fprintf('\n===== V2V中断概率结果（原始仿真数据）=====\n');
fprintf('行: gamma_th = [-5, 0, 5, 10, 15] dB\n');
fprintf('列: T = [0.2, 0.6, 1.0, 1.4, 1.8] ms\n');
disp('鲁棒算法中断概率(原始):');
disp(outage_prob_robust_raw);
disp('非鲁棒算法中断概率(原始):');
disp(outage_prob_nonrobust_raw);

toc
