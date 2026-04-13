% MAIN_V2I_CONVERGENCE_COMPLEXITY  脚本D：算法收敛性与计算复杂度分析
%
% 说明：
%   验证鲁棒算法（二分搜索）的收敛性，并对比鲁棒与非鲁棒算法的
%   计算复杂度（以迭代次数衡量）。
%
%   子图(a) — 收敛曲线：
%     横轴：二分搜索迭代次数（1, 2, 3, ...）
%     纵轴：V2I 容量（bps/Hz）
%     展示二分搜索如何逐步逼近最优功率值
%
%   子图(b) — 平均迭代次数随车辆密度变化：
%     横轴：车辆密度 N = 10, 20, 30, 40
%     纵轴：平均迭代次数
%     鲁棒算法约 28 次/配对，非鲁棒算法恒为 1 次/配对（闭式解）
%
% 输出图表：
%   V2I_Convergence_Curve.png / V2I_Convergence_Curve.pdf  — 收敛曲线
%   V2I_Avg_Iterations.png / V2I_Avg_Iterations.pdf         — 迭代次数对比
%
% 仿真模式：
%   fastMode = true：channNum_conv=50, channNum_time=30（约 8 秒）
%   fastMode = false：channNum_conv=200, channNum_time=50（约 3~5 分钟）
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
% fastMode = true：快速验证（~8秒，channNum_conv=50, channNum_time=30）
% fastMode = false：正式论文结果（~3-5分钟，channNum_conv=200, channNum_time=50）
fastMode = false;
if fastMode
    channNum_conv = 50;      % 收敛曲线的信道实现数
    channNum_time = 30;      % 复杂度统计的信道实现数
    fprintf('*** 快速测试模式 ***\n');
else
    channNum_conv = 200;    % 收敛曲线：充分多信道平均
    channNum_time = 50;      % 复杂度：足够测时间
    fprintf('*** 正式仿真模式 ***\n');
end
rng(3);  % 固定随机种子

%% =================================================================
%%  系统参数配置（与 main_V2V_Outage 保持一致）
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

r0 = 0.5;
dB_gamma0 = 5;
p0 = 1e-6;  % 与 sim_01/sim_03 保持一致

sig2 = 10^(dB_sig2 / 10);
gamma0 = 10^(dB_gamma0 / 10);
Pd_max = 10^(dB_Pd_max / 10);
Pc_max = 10^(dB_Pc_max / 10);

v_fixed = 60;               % 固定车速 60 km/h
T_fixed = 1;                 % 固定反馈周期 1 ms
N_list = [10, 20, 30, 40];  % 车辆密度列表（辆）

epsi_k = besselj(0, 2 * pi * (T_fixed * 1e-3) * (fc * 1e9) * (v_fixed / 3.6) / (3e8));
epsi_mk = epsi_k;

maxIter_limit = 50;  % 二分搜索最大迭代次数上限（安全保护）
d0 = sqrt(radius^2 - disBstoHwy^2);

%% =================================================================
%%  辅助函数：追踪鲁棒算法的收敛过程
%% =================================================================
% 本函数复制了 calOptPower.m 的核心逻辑，但额外记录每次迭代的
% V2I 容量和实际使用的 Pd 值，用于绘制收敛曲线。
%
% 与 calOptPower 的区别：
%   - calOptPower：只输出最终 (Pc_opt, Pd_opt)
%   - trackRobustConvergence：输出 (Pc_opt, Pd_opt) + 每次迭代的中间结果
%                             以及对应的 V2I 容量历史
%
% 输入：
%   g_mB, g_kB - V2I 链路和 V2V 干扰链路的信道增益（用于计算 V2I 容量）
%   maxIter    - 最大迭代次数上限
%
% 输出：
%   obj_history  - 列向量，每次迭代对应的 V2I 容量
%   Pd_history   - 列向量，每次迭代使用的 Pd 值
%   actual_iters - 标量，实际迭代次数

function [obj_history, Pd_history, actual_iters] = trackRobustConvergence( ...
    epsi, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0, ...
    g_mB, g_kB, maxIter)

    % -----------------------------------------------------------------
    % Step 1：计算闭式参考值 (Pc0, Pd0)
    %         用于判断功率空间落入 Case I / II / III 中的哪一个
    % -----------------------------------------------------------------
    den0 = alpha_mk * (1 - epsi_mk^2) * (1/p0 - 1) * epsi_k^2 * abs(h_k)^2 ...
           - (1 - epsi_k^2) * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
    Pc0 = (1 - epsi_k^2) * sig2 / den0;
    Pd0 = Pc0 * gamma0 * alpha_mk * (1 - epsi_mk^2) * (1 - p0) ...
          / (alpha_k * (1 - epsi_k^2) * p0);

    % 初始化历史记录向量（预分配 maxIter 大小的向量）
    obj_history = zeros(maxIter, 1);
    Pd_history = zeros(maxIter, 1);
    actual_iters = 0;

    % =================================================================
    %%  Case I：Pd_max <= Pd0（V2V 功率受限但可行区域充足）
    % =================================================================
    if Pd_max <= Pd0

        B = Pd_max * alpha_k * (1 - epsi_k^2);
        C = sig2 + Pc_max * epsi_mk^2 * alpha_mk * abs(h_mk)^2;
        D = Pc_max * alpha_mk * (1 - epsi_mk^2);
        tmp = 1 / (1 - p0) * exp(epsi_k^2 * abs(h_k)^2 / (1 - epsi_k^2));

        if exp(C * gamma0 / B) * (1 + D / B * gamma0) - tmp > 0
            % 分支 A：搜索最小 Pc（Pc_min），使约束刚好满足
            Pd_opt = Pd_max;
            P_left = 0; P_right = Pc_max;
            B = Pd_max * alpha_k * (1 - epsi_k^2);

            while abs(P_right - P_left) > epsi && actual_iters < maxIter
                actual_iters = actual_iters + 1;
                P_mid = (P_left + P_right) / 2;
                C = sig2 + P_mid * epsi_mk^2 * alpha_mk * abs(h_mk)^2;
                D = P_mid * alpha_mk * (1 - epsi_mk^2);

                % 记录本次迭代的 V2I 容量（基于当前 Pc = P_mid）
                obj_history(actual_iters) = log2(1 + P_mid * g_mB / (sig2 + Pd_opt * g_kB));
                Pd_history(actual_iters) = Pd_opt;

                if exp(C * gamma0 / B) * (1 + D / B * gamma0) - tmp > 0
                    P_right = P_mid;   % 约束可满足，Pc 可继续减小
                else
                    P_left = P_mid;   % 约束不满足，Pc 需要增大
                end
            end

        else
            % 分支 B：搜索最小 Pd（Pd_min），使约束刚好满足
            Pc_opt = Pc_max;
            P_left = 0; P_right = Pd_max;
            C = sig2 + Pc_max * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
            D = Pc_max * alpha_mk * (1 - epsi_mk^2);

            while abs(P_right - P_left) > epsi && actual_iters < maxIter
                actual_iters = actual_iters + 1;
                P_mid = (P_left + P_right) / 2;
                B = P_mid * alpha_k * (1 - epsi_k^2);

                % 记录本次迭代的 V2I 容量（基于当前 Pd = P_mid）
                obj_history(actual_iters) = log2(1 + Pc_opt * g_mB / (sig2 + P_mid * g_kB));
                Pd_history(actual_iters) = P_mid;

                if exp(C * gamma0 / B) * (1 + D / B * gamma0) - tmp < 0
                    P_right = P_mid;   % LHS 仍小，需要更大 Pd
                else
                    P_left = P_mid;   % LHS 已够大，Pd 可继续减小
                end
            end
        end

    % =================================================================
    %%  Case II：Pd_max > Pd0 但 Pc_max > Pc0（V2I 功率受限但可行）
    % =================================================================
    elseif Pc_max > Pc0

        num = (epsi_mk^2 * abs(h_mk)^2) / (1 - epsi_mk^2);
        A = Pd_max * alpha_k * epsi_k^2 * abs(h_k)^2;
        B = Pd_max * alpha_k * (1 - epsi_k^2);
        D = Pc_max * alpha_mk * (1 - epsi_mk^2);
        den1 = log(1 + B / (gamma0 * D));
        den2 = (A - sig2 * gamma0) / (gamma0 * D);

        if num - (den1 + den2) - log(p0) > 0
            % 分支 A：搜索最小 Pc
            Pd_opt = Pd_max;
            P_left = 0; P_right = Pc_max;
            A = Pd_max * alpha_k * epsi_k^2 * abs(h_k)^2;
            B = Pd_max * alpha_k * (1 - epsi_k^2);

            while abs(P_right - P_left) > epsi && actual_iters < maxIter
                actual_iters = actual_iters + 1;
                P_mid = (P_left + P_right) / 2;
                D = P_mid * alpha_mk * (1 - epsi_mk^2);
                den1 = log(1 + B / (gamma0 * D));
                den2 = (A - sig2 * gamma0) / (gamma0 * D);

                obj_history(actual_iters) = log2(1 + P_mid * g_mB / (sig2 + Pd_opt * g_kB));
                Pd_history(actual_iters) = Pd_opt;

                if num - (den1 + den2) - log(p0) > 0
                    P_right = P_mid;
                else
                    P_left = P_mid;
                end
            end

        else
            % 分支 B：搜索最小 Pd
            Pc_opt = Pc_max;
            P_left = 0; P_right = Pd_max;
            D = Pc_max * alpha_mk * (1 - epsi_mk^2);

            while abs(P_right - P_left) > epsi && actual_iters < maxIter
                actual_iters = actual_iters + 1;
                P_mid = (P_left + P_right) / 2;
                A = P_mid * alpha_k * epsi_k^2 * abs(h_k)^2;
                B = P_mid * alpha_k * (1 - epsi_k^2);

                den1 = log(1 + B / (gamma0 * D));
                den2 = (A - sig2 * gamma0) / (gamma0 * D);

                obj_history(actual_iters) = log2(1 + Pc_opt * g_mB / (sig2 + P_mid * g_kB));
                Pd_history(actual_iters) = P_mid;

                if num - (den1 + den2) - log(p0) < 0
                    P_right = P_mid;
                else
                    P_left = P_mid;
                end
            end
        end

    % =================================================================
    %%  Case III：不可行区域（Pd_max <= Pd0 且 Pc_max <= Pc0）
    %%             即使全功率仍无法满足 V2V 中断约束
    %%             回退：Pc=Pc_max，搜索最小 Pd
    % =================================================================
    else

        tmp = 1 / (1 - p0) * exp(epsi_k^2 * abs(h_k)^2 / (1 - epsi_k^2));
        Pc_opt = Pc_max;
        P_left = 0; P_right = Pd_max;
        C = sig2 + Pc_max * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
        D = Pc_max * alpha_mk * (1 - epsi_mk^2);

        while abs(P_right - P_left) > epsi && actual_iters < maxIter
            actual_iters = actual_iters + 1;
            P_mid = (P_left + P_right) / 2;
            B = P_mid * alpha_k * (1 - epsi_k^2);

            obj_history(actual_iters) = log2(1 + Pc_opt * g_mB / (sig2 + P_mid * g_kB));
            Pd_history(actual_iters) = P_mid;

            if exp(C * gamma0 / B) * (1 + D / B * gamma0) - tmp < 0
                P_right = P_mid;
            else
                P_left = P_mid;
            end
        end
    end

    % 截取实际使用的历史记录（去掉预分配的空余部分）
    if actual_iters < maxIter
        obj_history = obj_history(1 : actual_iters);
        Pd_history = Pd_history(1 : actual_iters);
    end

end


%% =================================================================
%%  子图1：收敛曲线
%% =================================================================
fprintf('正在生成收敛曲线...\n');

% 收集所有信道实现中各配对的收敛历史，然后做平均
all_obj_hist = {};
all_iters = zeros(channNum_conv, 1);       % 每信道实现的迭代次数
all_obj_nr = zeros(channNum_conv, 1);     % 每信道实现的非鲁棒 V2I 容量
all_obj_final_r = zeros(channNum_conv, 1);% 每信道实现的鲁棒收敛后 V2I 容量

for ch = 1 : channNum_conv

    % 生成车辆拓扑和信道（固定 numCUE=20, numDUE=20）
    [genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
        genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, ...
                     2.5 * v_fixed / 3.6, 20, 20);
    if genFlag == 1, continue; end

    % 大尺度衰落
    alpha_mB_c = zeros(1, 20);
    alpha_k_c = zeros(1, 20);
    alpha_kB_c = zeros(1, 20);
    alpha_mk_c = zeros(20, 20);

    for m = 1 : 20
        dist_mB = sqrt(vehPos(indCUE(m), 1)^2 + vehPos(indCUE(m), 2)^2);
        dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, fc) ...
                      + vehAntGain + bsAntGain - bsNoiseFigure;
        alpha_mB_c(m) = 10^(dB_alpha_mB / 10);

        for k = 1 : 20
            dist_mk = sqrt((vehPos(indCUE(m), 1) - vehPos(indDUE2(k), 1))^2 ...
                        + (vehPos(indCUE(m), 2) - vehPos(indDUE2(k), 2))^2);
            dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) ...
                          + 2 * vehAntGain - vehNoiseFigure;
            alpha_mk_c(m, k) = 10^(dB_alpha_mk / 10);
        end
    end

    for k = 1 : 20
        dist_k = sqrt((vehPos(indDUE(k), 1) - vehPos(indDUE2(k), 1))^2 ...
                    + (vehPos(indDUE(k), 2) - vehPos(indDUE2(k), 2))^2);
        dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) ...
                      + 2 * vehAntGain - vehNoiseFigure;
        alpha_k_c(k) = 10^(dB_alpha_k / 10);

        dist_kB = sqrt(vehPos(indDUE(k), 1)^2 + vehPos(indDUE(k), 2)^2);
        dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, fc) ...
                      + vehAntGain + bsAntGain - bsNoiseFigure;
        alpha_kB_c(k) = 10^(dB_alpha_kB / 10);
    end

    % 小尺度衰落
    h_mB_c = (randn(20, 1) + 1j * randn(20, 1)) / sqrt(2);
    h_k_c = (randn(20, 20) + 1j * randn(20, 20)) / sqrt(2);
    h_kB_c = (randn(20, 20) + 1j * randn(20, 20)) / sqrt(2);
    h_mk_c = (randn(20, 20) + 1j * randn(20, 20)) / sqrt(2);

    % 稀疏采样配对以减少计算量（每隔 5 个取一个）
    for m = 1 : 5 : 20
        for k = 1 : 5 : 20
            g_mB = alpha_mB_c(m) * abs(h_mB_c(m))^2;
            g_kB = alpha_kB_c(k) * abs(h_kB_c(m, k))^2;
            g_k = alpha_k_c(k) * abs(h_k_c(m, k))^2;
            if g_k <= 0, continue; end

            % 获取鲁棒算法的收敛历史（记录每次迭代的 V2I 容量）
            [obj_h, ~, iters] = trackRobustConvergence(1e-6, sig2, Pc_max, Pd_max, ...
                alpha_k_c(k), alpha_mk_c(m, k), epsi_k, epsi_mk, ...
                h_k_c(m, k), h_mk_c(m, k), p0, gamma0, g_mB, g_kB, maxIter_limit);

            all_obj_hist{end + 1} = obj_h;
            all_iters(length(all_obj_hist)) = iters;

            % 获取非鲁棒算法的 V2I 容量（闭式，一次计算完成）
            [Pd_nr, Pc_nr] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
                alpha_k_c(k), alpha_mk_c(m, k), h_k_c(m, k), h_mk_c(m, k), gamma0);
            all_obj_nr(length(all_obj_hist)) = log2(1 + Pc_nr * g_mB / (sig2 + Pd_nr * g_kB));
            all_obj_final_r(length(all_obj_hist)) = obj_h(iters);
        end
    end

    if mod(ch, 50) == 0
        fprintf('  收敛曲线进度: %d/%d\n', ch, channNum_conv);
    end
end

% 将不同长度的收敛历史对齐到相同长度，做平均
max_iters = max(all_iters);
mean_obj_hist = zeros(max_iters, 1);
count_hist = zeros(max_iters, 1);

for i = 1 : length(all_obj_hist)
    h = all_obj_hist{i};
    for j = 1 : length(h)
        mean_obj_hist(j) = mean_obj_hist(j) + h(j);
        count_hist(j) = count_hist(j) + 1;
    end
end
mean_obj_hist = mean_obj_hist ./ max(count_hist, 1);

mean_obj_nr = mean(all_obj_nr(all_obj_nr > 0));
mean_obj_final_r = mean(all_obj_final_r(all_obj_final_r > 0));
avg_iters = round(mean(all_iters(all_iters > 0)));

fprintf('收敛曲线: 平均迭代%.1f次, 鲁棒终值%.4f, 非鲁棒值%.4f\n', ...
    avg_iters, mean_obj_final_r, mean_obj_nr);


%% =================================================================
%%  子图2：平均迭代次数随车辆密度变化
%% =================================================================
fprintf('正在生成平均迭代次数数据...\n');

avg_iter_robust = zeros(1, length(N_list));
avg_iter_nonrobust = zeros(1, length(N_list));  % 恒为 1（闭式解）

for N_idx = 1 : length(N_list)
    N = N_list(N_idx);
    total_iter_r = 0;    % 鲁棒算法：总迭代次数累加
    pair_count = 0;     % 有效配对计数

    for ch = 1 : channNum_time
        [genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
            genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, ...
                         2.5 * v_fixed / 3.6, N, N);
        if genFlag == 1, continue; end

        % 大尺度衰落
        alpha_mB_N = zeros(1, N);
        alpha_k_N = zeros(1, N);
        alpha_kB_N = zeros(1, N);
        alpha_mk_N = zeros(N, N);

        for m = 1 : N
            dist_mB = sqrt(vehPos(indCUE(m), 1)^2 + vehPos(indCUE(m), 2)^2);
            dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, fc) ...
                          + vehAntGain + bsAntGain - bsNoiseFigure;
            alpha_mB_N(m) = 10^(dB_alpha_mB / 10);

            for k = 1 : N
                dist_mk = sqrt((vehPos(indCUE(m), 1) - vehPos(indDUE2(k), 1))^2 ...
                           + (vehPos(indCUE(m), 2) - vehPos(indDUE2(k), 2))^2);
                dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) ...
                              + 2 * vehAntGain - vehNoiseFigure;
                alpha_mk_N(m, k) = 10^(dB_alpha_mk / 10);
            end
        end

        for k = 1 : N
            dist_k = sqrt((vehPos(indDUE(k), 1) - vehPos(indDUE2(k), 1))^2 ...
                        + (vehPos(indDUE(k), 2) - vehPos(indDUE2(k), 2))^2);
            dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) ...
                         + 2 * vehAntGain - vehNoiseFigure;
            alpha_k_N(k) = 10^(dB_alpha_k / 10);

            dist_kB = sqrt(vehPos(indDUE(k), 1)^2 + vehPos(indDUE(k), 2)^2);
            dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, fc) ...
                          + vehAntGain + bsAntGain - bsNoiseFigure;
            alpha_kB_N(k) = 10^(dB_alpha_kB / 10);
        end

        % 小尺度衰落
        h_mB_N = (randn(N, 1) + 1j * randn(N, 1)) / sqrt(2);
        h_k_N = (randn(N, N) + 1j * randn(N, N)) / sqrt(2);
        h_kB_N = (randn(N, N) + 1j * randn(N, N)) / sqrt(2);
        h_mk_N = (randn(N, N) + 1j * randn(N, N)) / sqrt(2);

        % 遍历所有配对，统计迭代次数
        for m = 1 : N
            for k = 1 : N
                g_k_val = alpha_k_N(k) * abs(h_k_N(m, k))^2;
                if g_k_val <= 0, continue; end

                [~, ~, iters] = trackRobustConvergence(1e-6, sig2, Pc_max, Pd_max, ...
                    alpha_k_N(k), alpha_mk_N(m, k), epsi_k, epsi_mk, ...
                    h_k_N(m, k), h_mk_N(m, k), p0, gamma0, ...
                    alpha_mB_N(m) * abs(h_mB_N(m))^2, ...
                    alpha_kB_N(k) * abs(h_kB_N(m, k))^2, ...
                    maxIter_limit);

                total_iter_r = total_iter_r + iters;
                pair_count = pair_count + 1;
            end
        end
    end

    avg_iter_robust(N_idx) = total_iter_r / max(pair_count, 1);
    avg_iter_nonrobust(N_idx) = 1.0;  % 闭式解，恒为 1 次迭代
    fprintf('  N=%d: 鲁棒平均迭代%.1f次, 非鲁棒1次\n', N, avg_iter_robust(N_idx));
end


%% =================================================================
%%  绘图（学术规范，中文标注）
%% =================================================================
LineWidth = 1.5;
MarkerSize = 9;
FontSize = 12;
FontName = 'SimHei';  % 黑体以支持中文

fig = figure('Position', [100, 80, 1800, 650]);
set(gcf, 'Color', 'white', 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [32 12], 'PaperPosition', [0.5 0.5 31 11]);

% ==================================================================
%  子图1：收敛曲线
% ==================================================================
ax1 = subplot(1, 2, 1, 'Parent', fig);
set(ax1, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax1, 'on'); grid(ax1, 'minor'); hold(ax1, 'on');

% 绘制收敛曲线（鲁棒算法 V2I 容量随迭代次数的变化，逐次逼近最优值）
iter_axis = 1 : max_iters;
plot(ax1, iter_axis, mean_obj_hist, '-', ...
    'Color', [0.0 0.45 0.75], 'LineWidth', LineWidth + 0.5);

% 非鲁棒算法参考线（恒定值，闭式解无需迭代）
plot(ax1, [1, max_iters], [mean_obj_nr, mean_obj_nr], '--', ...
    'Color', [0.85 0.15 0.15], 'LineWidth', LineWidth);

% 鲁棒算法收敛终值参考线
plot(ax1, [1, max_iters], [mean_obj_final_r, mean_obj_final_r], ':', ...
    'Color', [0.0 0.6 0.3], 'LineWidth', LineWidth);

% 标注收敛点（平均迭代次数对应的 V2I 容量）
plot(ax1, avg_iters, mean_obj_final_r, '^k', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.0 0.45 0.75], 'LineWidth', 1.5);

xlabel('迭代次数', 'FontName', FontName, 'FontSize', FontSize + 1);
ylabel('V2I容量 (bps/Hz)', 'FontName', FontName, 'FontSize', FontSize + 1);
xlim([0, max_iters + 2]);
ylim([0, max(mean_obj_nr, mean_obj_final_r) * 1.25]);

legend(ax1, {'鲁棒算法（二分搜索）', '非鲁棒算法（闭式）', ...
    sprintf('鲁棒收敛值 (iter=%d)', avg_iters)}, ...
    'FontName', FontName, 'FontSize', FontSize - 1, ...
    'Location', 'southeast', 'Box', 'off');
title(ax1, '(a) 收敛曲线', 'FontName', FontName, 'FontSize', FontSize + 1);

% ==================================================================
%  子图2：平均迭代次数随车辆密度变化
% ==================================================================
ax2 = subplot(1, 2, 2, 'Parent', fig);
set(ax2, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax2, 'on'); grid(ax2, 'minor'); hold(ax2, 'on');

% 鲁棒算法：三角实线（迭代次数稳定，与密度弱相关）
plot(ax2, N_list, avg_iter_robust, '^-', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 3, 'MarkerFaceColor', [0.0 0.45 0.75], ...
    'Color', [0.0 0.45 0.75]);

% 非鲁棒算法：菱形虚线（恒为 1 次，闭式解）
plot(ax2, N_list, avg_iter_nonrobust, 'd--', 'LineWidth', LineWidth, ...
    'MarkerSize', MarkerSize + 3, 'MarkerFaceColor', [0.85 0.15 0.15], ...
    'Color', [0.85 0.15 0.15]);
% 为非鲁棒算法添加闭式解标注
text(N_list(end) + 2, avg_iter_nonrobust(end), '闭式解（无迭代）', ...
    'FontName', FontName, 'FontSize', FontSize - 1, ...
    'Color', [0.85 0.15 0.15], 'VerticalAlignment', 'middle');

xlabel('车辆密度 N', 'FontName', FontName, 'FontSize', FontSize + 1);
ylabel('平均迭代次数', 'FontName', FontName, 'FontSize', FontSize + 1);
xlim([5, 50]);
ylim([-2, max(avg_iter_robust) * 1.3]);  % Y轴从-2开始，凸显非鲁棒算法的零值

legend(ax2, {'鲁棒算法', '非鲁棒算法'}, ...
    'FontName', FontName, 'FontSize', FontSize - 1, ...
    'Location', 'northwest', 'Box', 'off');
title(ax2, {'(b) 平均迭代次数随密度变化'}, 'FontName', FontName, 'FontSize', FontSize + 1);

% 输出图像
print('-dpng', '-r300', fullfile(outputFolder, 'sim_04_Algorithm_Convergence_Complexity.png'));
print('-dpdf', '-r300', fullfile(outputFolder, 'sim_04_Algorithm_Convergence_Complexity.pdf'));
fprintf('图已保存: %s/sim_04_Algorithm_Convergence_Complexity.png / .pdf\n', outputFolder);

fprintf('\n=== 仿真完成 ===\n');
toc
