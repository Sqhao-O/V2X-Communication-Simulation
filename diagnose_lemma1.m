% DIAGNOSE_LEMMA1  测试 Lemma 1 精确中断概率约束的效果
%
% 目的：
%   验证论文中 Lemma 1 给出的闭式中断概率上界公式的准确性。
%   通过蒙特卡洛仿真对比：
%     - 鲁棒算法：基于 Lemma 1 的 Markov 上界设计功率
%     - 非鲁棒算法：忽略信道误差，基于估计信道设计功率
%     在同一个随机信道实现下的实际中断概率
%
% Lemma 1 回顾（论文公式 10/11）：
%   给定功率分配 (Pc, Pd) 和信道估计 (h_k, h_mk)，
%   实际 SINR 低于阈值 gamma0 的概率具有以下闭式上界：
%
%   定义中间变量：
%     A = Pd * alpha_k * epsi_k^2 * |h_k|^2      （V2V 直连确定分量）
%     B = Pd * alpha_k * (1 - epsi_k^2)           （V2V 直连随机分量强度）
%     C = sig2 + Pc * alpha_mk * epsi_mk^2 * |h_mk|^2  （干扰确定分量）
%     D = Pc * alpha_mk * (1 - epsi_mk^2)         （干扰随机分量强度）
%
%   Case I (C*gamma0 >= A)：Pr(SINR < gamma0) = 1 - exp(-(C*gamma0-A)/B)/(1+D*gamma0/B)
%   Case II (其他)：         Pr(SINR < gamma0) = exp((A-C*gamma0)/(gamma0*D))/(1+B/(gamma0*D))
%
% 仿真逻辑：
%   1. 对给定信道估计 (h_k, h_mk)，计算鲁棒和非鲁棒功率分配
%   2. 生成大量随机信道误差样本（模拟延迟 CSI 的不确定性）
%   3. 统计实际 SINR < gamma0 的频率 → 实际中断概率
%   4. 与 Lemma 1 的理论预测值对比
%
% 预期结果：
%   - 鲁棒算法的实际中断概率 ≤ p0（约束有效）
%   - 非鲁棒算法的实际中断概率 >> p0（严重超限）
%   - Lemma 1 的理论值与蒙特卡洛结果高度一致
%
% 被测函数：calOptPower.m, calOptPower_nonrobust.m
% -----------------------------------------------------------------

clear; clc;

%% =================================================================
%%  系统参数设置
%% =================================================================
sig2 = 10^(-114/10);     % 噪声功率（线性值，W）
gamma0 = 10^(5/10);      % V2V SINR 阈值（5 dB）
Pd_max = 10^(23/10);     % V2V 最大功率（23 dBm，线性值）
Pc_max = 10^(23/10);     % V2I 最大功率（23 dBm，线性值）

% 测试不同的 p0 值，观察鲁棒算法对约束收紧的响应
p0_list = [1e-3, 1e-4, 1e-5];

% 信道时间相关系数（v=60km/h, T=1ms, fc=2GHz）
epsi_k = besselj(0, 2 * pi * (1e-3) * (2e9) * (60/3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('gamma0 = %.2f, epsi_k = %.4f\n', gamma0, epsi_k);
fprintf('Pd_max = %.2f W, Pc_max = %.2f W\n', Pd_max, Pc_max);

%% =================================================================
%%  信道条件：典型 V2V 场景
%% =================================================================
% alpha_k ≈ 10^(-6)：V2V 直连链路约 20m 距离，含路径损耗 + 阴影 + 天线增益
% alpha_mk ≈ 10^(-8)：V2I→V2V 干扰链路，CUE 距 DUE 接收端较远（约 100m）
alpha_k = 1e-6;
alpha_mk = 1e-8;

fprintf('\n=== 信道状态 ===\n');
fprintf('alpha_k = %.2e, alpha_mk = %.2e\n', alpha_k, alpha_mk);

%% =================================================================
%%  诊断 1：对不同 p0 值进行蒙特卡洛中断概率测试
%% =================================================================
% 对每个 p0：
%   - numTest = 50 次独立信道估计采样
%   - 每次估计后用 numErr = 500 次误差采样评估实际中断概率
%   - 总样本数 = 50 × 500 = 25,000
for p0 = p0_list
    fprintf('\n========== p0 = %.0e ==========\n', p0);

    %% 蒙特卡洛测试
    numTest = 50;    % 信道估计采样数
    numErr = 500;    % 每次估计的误差采样数
    outage_count_robust = 0;
    outage_count_nonrobust = 0;
    total_samples = 0;

    for i = 1:numTest
        % 随机生成信道估计值（复高斯 CN(0,1)）
        h_k = (randn + 1j*randn) / sqrt(2);
        h_mk = (randn + 1j*randn) / sqrt(2);

        % 鲁棒功率分配：基于 Lemma 1 的 Markov 上界
        [Pd_robust, Pc_robust] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
            alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

        % 非鲁棒功率分配：假设完美 CSI (epsi = 1)
        [Pd_nonrobust, Pc_nonrobust] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
            alpha_k, alpha_mk, h_k, h_mk, gamma0);

        % 测试实际中断概率（对信道误差做蒙特卡洛）
        for j = 1:numErr
            % 生成信道估计误差（零均值复高斯，方差 = 1-epsi^2）
            e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
            e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);

            % 构造实际信道 = 估计值×相关系数 + 误差项
            hk_actual = epsi_k * h_k + e_k;
            hmk_actual = epsi_mk * h_mk + e_mk;
            gk_actual = alpha_k * abs(hk_actual)^2;
            gmk_actual = alpha_mk * abs(hmk_actual)^2;

            % 计算实际 SINR
            sinr_robust = Pd_robust * gk_actual / (sig2 + Pc_robust * gmk_actual);
            sinr_nonrobust = Pd_nonrobust * gk_actual / (sig2 + Pc_nonrobust * gmk_actual);

            % 中断判断：SINR < gamma0
            if sinr_robust < gamma0
                outage_count_robust = outage_count_robust + 1;
            end
            if sinr_nonrobust < gamma0
                outage_count_nonrobust = outage_count_nonrobust + 1;
            end
            total_samples = total_samples + 1;
        end
    end

    fprintf('鲁棒算法实际中断概率:   %.4f (%.2f%%)\n', outage_count_robust/total_samples, 100*outage_count_robust/total_samples);
    fprintf('非鲁棒算法实际中断概率: %.4f (%.2f%%)\n', outage_count_nonrobust/total_samples, 100*outage_count_nonrobust/total_samples);
    fprintf('目标p0:                 %.4f (%.2f%%)\n', p0, 100*p0);
end

%% =================================================================
%%  诊断 2：单个案例的详细分析（p0 = 1e-4）
%% =================================================================
% 选取一个具体信道实现，计算 Lemma 1 的各个参数，
% 验证理论公式的数值计算过程
fprintf('\n========== 单个案例分析 (p0=1e-4) ==========\n');
p0 = 1e-4;
h_k = (randn + 1j*randn) / sqrt(2);
h_mk = (randn + 1j*randn) / sqrt(2);

[Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

fprintf('分配的功率:\n');
fprintf('Pd_opt = %.4f W (%.1f%% of max)\n', Pd_opt, 100*Pd_opt/Pd_max);
fprintf('Pc_opt = %.4f W (%.1f%% of max)\n', Pc_opt, 100*Pc_opt/Pc_max);

% ---- 计算 Lemma 1 所需的中间变量 ----
% A: V2V 直连链路的确定分量（来自信道估计的可预测部分）
A = Pd_opt * alpha_k * epsi_k^2 * abs(h_k)^2;
% B: V2V 直连链路的随机分量强度（来自信道估计误差）
B = Pd_opt * alpha_k * (1 - epsi_k^2);
% C: 干扰 + 噪声的确定分量
C = sig2 + Pc_opt * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
% D: 干扰的随机分量强度
D = Pc_opt * alpha_mk * (1 - epsi_mk^2);

fprintf('\nLemma 1参数:\n');
fprintf('A = %.6e (V2V直连确定分量)\n', A);
fprintf('B = %.6e (V2V直连随机分量)\n', B);
fprintf('C = %.6e (干扰+噪声确定分量)\n', C);
fprintf('D = %.6e (干扰随机分量)\n', D);
fprintf('C*gamma0 = %.6e, A = %.6e\n', C*gamma0, A);

% 根据 Case 类型选择对应的闭式公式
if C * gamma0 >= A
    % Case I：确定分量主导（干扰确定分量乘以阈值大于信号确定分量）
    % Pr(SINR < gamma0) = 1 - exp(-(C*gamma0 - A)/B) / (1 + D*gamma0/B)
    fprintf('Case I: 使用公式(10)\n');
    p_theory = 1 - exp(-(C*gamma0 - A)/B) / (1 + D*gamma0/B);
else
    % Case II：随机分量主导
    % Pr(SINR < gamma0) = exp((A - C*gamma0)/(gamma0*D)) / (1 + B/(gamma0*D))
    fprintf('Case II: 使用公式(11)\n');
    p_theory = exp((A - C*gamma0)/(gamma0*D)) / (1 + B/(gamma0*D));
end
fprintf('理论中断概率 = %.6e (%.4f%%)\n', p_theory, 100*p_theory);
fprintf('目标p0 = %.6e (%.4f%%)\n', p0, 100*p0);
fprintf('理论值 vs 目标: 理论值应 <= 目标值（鲁棒约束有效）\n');
