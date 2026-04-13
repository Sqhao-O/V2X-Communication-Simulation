% DIAGNOSE_LEMMA1  测试Lemma 1精确中断概率约束的效果

clear; clc;

%% 参数设置
sig2 = 10^(-114/10);
gamma0 = 10^(5/10);
Pd_max = 10^(23/10);
Pc_max = 10^(23/10);

% 测试不同的p0值
p0_list = [1e-3, 1e-4, 1e-5];

epsi_k = besselj(0, 2 * pi * (1e-3) * (2e9) * (60/3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('gamma0 = %.2f, epsi_k = %.4f\n', gamma0, epsi_k);
fprintf('Pd_max = %.2f W, Pc_max = %.2f W\n', Pd_max, Pc_max);

%% 测试信道条件
alpha_k = 1e-6;
alpha_mk = 1e-8;

fprintf('\n=== 信道状态 ===\n');
fprintf('alpha_k = %.2e, alpha_mk = %.2e\n', alpha_k, alpha_mk);

%% 对每个p0值进行测试
for p0 = p0_list
    fprintf('\n========== p0 = %.0e ==========\n', p0);

    %% 蒙特卡洛测试
    numTest = 50;
    numErr = 500;
    outage_count_robust = 0;
    outage_count_nonrobust = 0;
    total_samples = 0;

    for i = 1:numTest
        h_k = (randn + 1j*randn) / sqrt(2);
        h_mk = (randn + 1j*randn) / sqrt(2);

        % 鲁棒功率分配
        [Pd_robust, Pc_robust] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
            alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

        % 非鲁棒功率分配
        [Pd_nonrobust, Pc_nonrobust] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
            alpha_k, alpha_mk, h_k, h_mk, gamma0);

        % 测试中断概率
        for j = 1:numErr
            e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
            e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);

            hk_actual = epsi_k * h_k + e_k;
            hmk_actual = epsi_mk * h_mk + e_mk;
            gk_actual = alpha_k * abs(hk_actual)^2;
            gmk_actual = alpha_mk * abs(hmk_actual)^2;

            sinr_robust = Pd_robust * gk_actual / (sig2 + Pc_robust * gmk_actual);
            sinr_nonrobust = Pd_nonrobust * gk_actual / (sig2 + Pc_nonrobust * gmk_actual);

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

%% 详细分析一个案例
fprintf('\n========== 单个案例分析 (p0=1e-4) ==========\n');
p0 = 1e-4;
h_k = (randn + 1j*randn) / sqrt(2);
h_mk = (randn + 1j*randn) / sqrt(2);

[Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

fprintf('分配的功率:\n');
fprintf('Pd_opt = %.4f W (%.1f%% of max)\n', Pd_opt, 100*Pd_opt/Pd_max);
fprintf('Pc_opt = %.4f W (%.1f%% of max)\n', Pc_opt, 100*Pc_opt/Pc_max);

% 计算理论中断概率
A = Pd_opt * alpha_k * epsi_k^2 * abs(h_k)^2;
B = Pd_opt * alpha_k * (1 - epsi_k^2);
C = sig2 + Pc_opt * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
D = Pc_opt * alpha_mk * (1 - epsi_mk^2);

fprintf('\nLemma 1参数:\n');
fprintf('A = %.6e, B = %.6e\n', A, B);
fprintf('C = %.6e, D = %.6e\n', C, D);
fprintf('C*gamma0 = %.6e, A = %.6e\n', C*gamma0, A);

if C * gamma0 >= A
    fprintf('Case I: 使用公式(10)\n');
    p_theory = 1 - exp(-(C*gamma0 - A)/B) / (1 + D*gamma0/B);
else
    fprintf('Case II: 使用公式(11)\n');
    p_theory = exp((A - C*gamma0)/(gamma0*D)) / (1 + B/(gamma0*D));
end
fprintf('理论中断概率 = %.6e (%.4f%%)\n', p_theory, 100*p_theory);
fprintf('目标p0 = %.6e (%.4f%%)\n', p0, 100*p0);
