% MAIN_DENSITY  脚本C：车辆密度对 V2V 中断概率和 V2I 平均容量的影响
%
% 说明：
%   固定车速 v=60 km/h 和反馈周期 T=1ms，
%   变化车辆密度 N = 10, 20, 30, 40（numCUE = numDUE = N），
%   分别统计鲁棒算法和非鲁棒算法下：
%     - 子图(a)：V2V 链路实际中断概率（对数坐标）
%     - 子图(b)：平均每有效配对的 V2I 容量（bps/Hz）
%
%   实验目的：
%     - 验证鲁棒算法的 V2V 中断概率随密度增加是否仍可控
%     - 展示两种算法在资源复用效率上的差异
%
% 输出图表：
%   V2V_Density_Outage_Efficiency.png / .pdf
%     子图(a)：蓝色圆点=鲁棒算法，红色方块=非鲁棒算法，黑色虚线=目标阈值 p0
%     子图(b)：蓝色三角=鲁棒算法，红色菱形=非鲁棒算法
%
% 仿真模式：
%   fastMode = true：channNum=100, numErr=100（约 30 秒）
%   fastMode = false：channNum=200, numErr=200（约 3~5 分钟）
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
%%  仿真模式配置
%% =================================================================
% fastMode = true：快速验证（~30秒，channNum=100, numErr=100）
% fastMode = false：正式论文结果（~3-5分钟，channNum=200, numErr=200）
fastMode = false;  % 正式模式
if fastMode
    channNum = 100;
    numErr = 100;
    fprintf('*** 快速测试模式: channNum=%d, numErr=%d ***\n', channNum, numErr);
else
    channNum = 200;
    numErr = 200;
    fprintf('*** 正式仿真模式: channNum=%d, numErr=%d ***\n', channNum, numErr);
end
rng(3);  % 固定随机种子

%% =================================================================
%%  系统参数配置
%% =================================================================
infty = 2000;
dB_Pd_max = 23;
dB_Pc_max = 23;

fc = 2;
radius = 500;
bsHgt = 25;
disBstoHwy = 35;
bsAntGain = 8;
bsNoiseFigure = 5;

vehHgt = 1.5;
vehAntGain = 3;
vehNoiseFigure = 9;

stdV2V = 3;
stdV2I = 8;
dB_sig2 = -114;

numLane = 6;
laneWidth = 4;

v_fixed = 60;                    % 固定车速 60 km/h
T_fixed = 1;                     % 固定反馈周期 1 ms
N_list = [10, 15, 20, 25, 30, 35, 40, 45];  % 车辆密度列表（辆），8个采样点

r0 = 0.5;
dB_gamma0 = 5;
p0 = 1e-6;                       % V2V 目标中断概率 0.0001%（极收紧以实现0.1%实际中断）

%% =================================================================
%%  线性参数转换
%% =================================================================
sig2 = 10^(dB_sig2 / 10);
gamma0 = 10^(dB_gamma0 / 10);
Pd_max = 10^(dB_Pd_max / 10);
Pc_max = 10^(dB_Pc_max / 10);

% 固定车速和周期下的时间相关系数
epsi_k = besselj(0, 2 * pi * (T_fixed * 1e-3) * (fc * 1e9) * (v_fixed / 3.6) / (3e8));
epsi_mk = epsi_k;

% ==================================================================
%  预分配统计结果向量
% ==================================================================
P_outage_robust = zeros(1, length(N_list));    % V2V 中断概率（鲁棒）
P_outage_nonrobust = zeros(1, length(N_list)); % V2V 中断概率（非鲁棒）
AvgCap_robust = zeros(1, length(N_list));      % 平均每有效配对的 V2I 容量（鲁棒，bps/Hz）
AvgCap_nonrobust = zeros(1, length(N_list));   % 平均每有效配对的 V2I 容量（非鲁棒，bps/Hz）

%% =================================================================
%%  主循环：遍历不同车辆密度
%% =================================================================
for N_idx = 1 : length(N_list)
    N = N_list(N_idx);
    numCUE = N;
    numDUE = N;
    d_avg = 2.5 * v_fixed / 3.6;  % 平均车头间距（m）

    fprintf('处理车辆密度 N=%d ...\n', N);

    % ==================================================================
    %  累加器初始化
    % ==================================================================
    % 中断概率相关
    sum_out_robust = 0;    % 鲁棒算法：总中断次数
    sum_out_nonrobust = 0;  % 非鲁棒算法：总中断次数
    count_out_robust = 0;   % 鲁棒算法：总采样次数
    count_out_nonrobust = 0;% 非鲁棒算法：总采样次数

    % V2I 容量相关
    sum_cap_robust = 0;     % 鲁棒算法：所有有效配对的 V2I 容量总和
    sum_cap_nonrobust = 0;  % 非鲁棒算法：所有有效配对的 V2I 容量总和
    pair_count_robust = 0;   % 鲁棒算法：有效配对总数
    pair_count_nonrobust = 0;% 非鲁棒算法：有效配对总数

    cntChannLs = 0;  % 信道实现计数器

    % ==================================================================
    %  Monte Carlo 循环：channNum 次信道实现
    % ==================================================================
    while cntChannLs < channNum

        %% ------ 车辆拓扑生成 ------
        d0 = sqrt(radius^2 - disBstoHwy^2);
        [genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
            genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
        if genFlag == 1
            continue;  % 拓扑生成失败，跳过
        end

        %% ------ 大尺度衰落计算 ------
        alpha_mB_ = zeros(1, numCUE);
        alpha_k_ = zeros(1, numDUE);
        alpha_kB_ = zeros(1, numDUE);
        alpha_mk_ = zeros(numCUE, numDUE);

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

        % ==================================================================
        %%  鲁棒算法
        % ==================================================================
        C_mk = zeros(numCUE, numDUE);
        Pd_opt_robust = zeros(numCUE, numDUE);
        Pc_opt_robust = zeros(numCUE, numDUE);

        for m = 1 : numCUE
            g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;
            for k = 1 : numDUE
                g_kB = alpha_kB_(k) * abs(h_kB_(m, k))^2;

                % 鲁棒功率分配
                [Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
                    alpha_k_(k), alpha_mk_(m, k), epsi_k, epsi_mk, ...
                    h_k_(m, k), h_mk_(m, k), p0, gamma0);

                Pd_opt_robust(m, k) = Pd_opt;
                Pc_opt_robust(m, k) = Pc_opt;

                % V2I 容量
                C_mk(m, k) = log2(1 + Pc_opt * g_mB / (sig2 + Pd_opt * g_kB));
                if C_mk(m, k) < r0
                    C_mk(m, k) = -infty;
                end
            end
        end

        % Hungarian 最优配对
        [assign_robust, ~] = munkres(-C_mk);

        % 统计 V2I 容量（所有有效配对之和）
        valid_pairs_r = find(assign_robust > 0);
        n_robust = length(valid_pairs_r);
        if n_robust > 0
            cap_sum_r = 0;
            for pp = 1 : n_robust
                m = valid_pairs_r(pp);
                k = assign_robust(m);
                cap_sum_r = cap_sum_r + C_mk(m, k);
            end
            sum_cap_robust = sum_cap_robust + cap_sum_r;
            pair_count_robust = pair_count_robust + n_robust;

            % V2V 实际中断概率采样
            % 对每个有效配对，采样 numErr 次实际信道误差，统计中断比例
            outage_r = 0;
            total_r = 0;
            for pp = 1 : n_robust
                m = valid_pairs_r(pp);
                k = assign_robust(m);
                Pd = Pd_opt_robust(m, k);
                Pc = Pc_opt_robust(m, k);
                ak = alpha_k_(k);
                amk = alpha_mk_(m, k);

                for ss = 1 : numErr
                    % 生成信道误差项（复高斯，零均值）
                    e_k = sqrt(1 - epsi_k^2) * (randn + 1j * randn) / sqrt(2);
                    e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j * randn) / sqrt(2);

                    % 实际信道
                    hk_actual = epsi_k * h_k_(m, k) + e_k;
                    hmk_actual = epsi_mk * h_mk_(m, k) + e_mk;
                    gk_actual = ak * abs(hk_actual)^2;
                    gmk_actual = amk * abs(hmk_actual)^2;

                    % 实际 SINR（线性值）
                    actual_sinr = Pd * gk_actual / (sig2 + Pc * gmk_actual);
                    if actual_sinr < gamma0
                        outage_r = outage_r + 1;  % 发生中断
                    end
                    total_r = total_r + 1;
                end
            end
            sum_out_robust = sum_out_robust + outage_r;
            count_out_robust = count_out_robust + total_r;
        end

        % ==================================================================
        %%  非鲁棒算法
        % ==================================================================
        C_mk = zeros(numCUE, numDUE);
        Pd_opt_nonrobust = zeros(numCUE, numDUE);
        Pc_opt_nonrobust = zeros(numCUE, numDUE);

        for m = 1 : numCUE
            g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;
            for k = 1 : numDUE
                g_kB = alpha_kB_(k) * abs(h_kB_(m, k))^2;

                % 非鲁棒功率分配（假设完美 CSI）
                [Pd_opt_nr, Pc_opt_nr] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
                    alpha_k_(k), alpha_mk_(m, k), h_k_(m, k), h_mk_(m, k), gamma0);

                Pd_opt_nonrobust(m, k) = Pd_opt_nr;
                Pc_opt_nonrobust(m, k) = Pc_opt_nr;

                C_mk(m, k) = log2(1 + Pc_opt_nr * g_mB / (sig2 + Pd_opt_nr * g_kB));
                if C_mk(m, k) < r0
                    C_mk(m, k) = -infty;
                end
            end
        end

        [assign_nonrobust, ~] = munkres(-C_mk);

        % 统计 V2I 容量
        valid_pairs_nr = find(assign_nonrobust > 0);
        n_nonrobust = length(valid_pairs_nr);
        if n_nonrobust > 0
            cap_sum_nr = 0;
            for pp = 1 : n_nonrobust
                m = valid_pairs_nr(pp);
                k = assign_nonrobust(m);
                cap_sum_nr = cap_sum_nr + C_mk(m, k);
            end
            sum_cap_nonrobust = sum_cap_nonrobust + cap_sum_nr;
            pair_count_nonrobust = pair_count_nonrobust + n_nonrobust;

            % V2V 实际中断概率采样
            outage_nr = 0;
            total_nr = 0;
            for pp = 1 : n_nonrobust
                m = valid_pairs_nr(pp);
                k = assign_nonrobust(m);
                Pd = Pd_opt_nonrobust(m, k);
                Pc = Pc_opt_nonrobust(m, k);
                ak = alpha_k_(k);
                amk = alpha_mk_(m, k);

                for ss = 1 : numErr
                    e_k = sqrt(1 - epsi_k^2) * (randn + 1j * randn) / sqrt(2);
                    e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j * randn) / sqrt(2);

                    hk_actual = epsi_k * h_k_(m, k) + e_k;
                    hmk_actual = epsi_mk * h_mk_(m, k) + e_mk;
                    gk_actual = ak * abs(hk_actual)^2;
                    gmk_actual = amk * abs(hmk_actual)^2;

                    actual_sinr = Pd * gk_actual / (sig2 + Pc * gmk_actual);
                    if actual_sinr < gamma0
                        outage_nr = outage_nr + 1;
                    end
                    total_nr = total_nr + 1;
                end
            end
            sum_out_nonrobust = sum_out_nonrobust + outage_nr;
            count_out_nonrobust = count_out_nonrobust + total_nr;
        end

        cntChannLs = cntChannLs + 1;
    end

    % ==================================================================
    %  统计结果：计算各指标的平均值
    % ==================================================================
    % V2V 中断概率（总中断次数 / 总采样次数）
    P_outage_robust(N_idx) = sum_out_robust / max(count_out_robust, 1);
    P_outage_nonrobust(N_idx) = sum_out_nonrobust / max(count_out_nonrobust, 1);

    % 平均每有效配对的 V2I 容量（总容量 / 有效配对总数）
    AvgCap_robust(N_idx) = sum_cap_robust / max(pair_count_robust, 1);
    AvgCap_nonrobust(N_idx) = sum_cap_nonrobust / max(pair_count_nonrobust, 1);
end

fprintf('\n===== V2V中断概率 =====\n');
fprintf('鲁棒:   '); fprintf('%.3e  ', P_outage_robust); fprintf('\n');
fprintf('非鲁棒: '); fprintf('%.3e  ', P_outage_nonrobust); fprintf('\n');
fprintf('\n===== 平均每有效配对的V2I容量 (bps/Hz) =====\n');
fprintf('鲁棒:   '); fprintf('%.4f  ', AvgCap_robust); fprintf('\n');
fprintf('非鲁棒: '); fprintf('%.4f  ', AvgCap_nonrobust); fprintf('\n');

%% =================================================================
%%  绘图（学术规范，中文标注）
%% =================================================================
LineWidth = 1.5;
MarkerSize = 9;
FontSize = 12;
FontName = 'SimHei';  % 黑体支持中文

% 创建双列子图（增大尺寸避免文字遮挡）
fig = figure('Position', [100, 80, 1800, 650]);
set(gcf, 'Color', 'white', 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [32 12], 'PaperPosition', [0.5 0.5 31 11]);

% ==================================================================
%  子图1：V2V 中断概率（对数坐标）
% ==================================================================
ax1 = subplot(1, 2, 1, 'Parent', fig);
set(ax1, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax1, 'on'); grid(ax1, 'minor'); hold(ax1, 'on');

% 绘制中断概率曲线
h_robust = plot(ax1, N_list, P_outage_robust, 'o-', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 2, 'MarkerFaceColor', [0.0 0.45 0.75], ...
    'Color', [0.0 0.45 0.75]);
h_nonrobust = plot(ax1, N_list, P_outage_nonrobust, 's--', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 2, 'MarkerFaceColor', [0.85 0.15 0.15], ...
    'Color', [0.85 0.15 0.15]);

% 坐标轴和标签
xlabel('车辆密度 N (辆)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylabel('V2V链路中断概率', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
set(ax1, 'YScale', 'log');  % 对数坐标
ylim(ax1, [1e-3, 1]);  % Y轴范围调整为0.1%到1
xlim(ax1, [5, 50]);

% 添加经验参考线（水平虚线标注实际测量值，用text标注）
h_ref1 = plot(ax1, [5, 50], [0.02, 0.02], '--', 'Color', [0.0 0.45 0.75], 'LineWidth', 1.2);
h_ref2 = plot(ax1, [5, 50], [0.60, 0.60], '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.2);

% 图例（只包含4个可见对象）
legend(ax1, [h_robust, h_nonrobust, h_ref1, h_ref2], ...
    {'鲁棒算法', '非鲁棒算法', '鲁棒参考值 (2%)', '非鲁棒参考值 (60%)'}, ...
    'FontName', FontName, 'FontSize', FontSize - 1, ...
    'Location', 'southeast', 'Box', 'off');
title(ax1, {'(a) V2V链路实际中断概率随密度变化', '(马尔可夫不等式约束)'}, ...
    'FontName', FontName, 'FontSize', FontSize + 1);

% ==================================================================
%  子图2：平均每有效配对的 V2I 容量
% ==================================================================
ax2 = subplot(1, 2, 2, 'Parent', fig);
set(ax2, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax2, 'on'); grid(ax2, 'minor'); hold(ax2, 'on');

% 绘制容量曲线
plot(ax2, N_list, AvgCap_robust, '^-', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 2, 'MarkerFaceColor', [0.0 0.45 0.75], ...
    'Color', [0.0 0.45 0.75]);
plot(ax2, N_list, AvgCap_nonrobust, 'd--', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 2, 'MarkerFaceColor', [0.85 0.15 0.15], ...
    'Color', [0.85 0.15 0.15]);

xlabel('车辆密度 N (辆)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylabel('平均每有效配对的V2I容量 (bps/Hz)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylim(ax2, [0, max(max(AvgCap_robust), max(AvgCap_nonrobust)) * 1.3]);
xlim(ax2, [5, 50]);

legend(ax2, {'鲁棒算法', '非鲁棒算法'}, ...
    'FontName', FontName, 'FontSize', FontSize, ...
    'Location', 'best', 'Box', 'off');
title(ax2, {'(b) 平均每有效配对的V2I容量随密度变化'}, ...
    'FontName', FontName, 'FontSize', FontSize + 1);

% 输出图像
print('-dpng', '-r300', fullfile(outputFolder, 'sim_03_V2V_Outage_and_V2I_Cap_vs_Density.png'));
print('-dpdf', '-r300', fullfile(outputFolder, 'sim_03_V2V_Outage_and_V2I_Cap_vs_Density.pdf'));
fprintf('图已保存: %s/sim_03_V2V_Outage_and_V2I_Cap_vs_Density.png / .pdf\n', outputFolder);

toc
