% DIAGNOSE_OUTAGE  诊断V2V中断概率过高的原因
%
% 本脚本用于分析calOptPower函数的实际行为，检查约束是否被满足

clear; clc;

%% 参数设置（与sim_03一致）
dB_Pd_max = 23; dB_Pc_max = 23;
fc = 2; radius = 500; bsHgt = 25; disBstoHwy = 35;
bsAntGain = 8; bsNoiseFigure = 5;
vehHgt = 1.5; vehAntGain = 3; vehNoiseFigure = 9;
stdV2V = 3; stdV2I = 8;
dB_sig2 = -114;

v_fixed = 60;  % km/h
T_fixed = 1;   % ms

r0 = 0.5;
dB_gamma0 = 5;
p0 = 1e-15;  % 极大收紧以接近0.1%目标

%% 线性转换
sig2 = 10^(dB_sig2 / 10);
gamma0 = 10^(dB_gamma0 / 10);
Pd_max = 10^(dB_Pd_max / 10);
Pc_max = 10^(dB_Pc_max / 10);

% 时间相关系数
epsi_k = besselj(0, 2 * pi * (T_fixed * 1e-3) * (fc * 1e9) * (v_fixed / 3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('p0 = %.0e (目标中断概率)\n', p0);
fprintf('gamma0 = %.2f (%.1f dB, 线性SINR阈值)\n', gamma0, dB_gamma0);
fprintf('epsi_k = %.4f (时间相关系数)\n', epsi_k);
fprintf('Pd_max = %.2f W (%.0f dBm)\n', Pd_max, dB_Pd_max);
fprintf('Pc_max = %.2f W (%.0f dBm)\n', Pc_max, dB_Pc_max);
fprintf('sig2 = %.2e W (%.0f dBm)\n', sig2, dB_sig2);

%% 模拟一些典型信道条件
fprintf('\n=== 测试不同信道条件下的功率分配 ===\n');

% 典型大尺度衰落值（根据genPL计算，V2V距离约20-50m）
alpha_k_test = 1e-6;   % V2V直连链路 ~ -60 dB
alpha_mk_test = 1e-8;  % V2I->V2V干扰链路 ~ -80 dB (更远)

fprintf('alpha_k = %.2e (V2V直连)\n', alpha_k_test);
fprintf('alpha_mk = %.2e (V2I->V2V干扰)\n', alpha_mk_test);

% 测试多个随机信道实现
numTest = 100;
outage_count_robust = 0;
outage_count_nonrobust = 0;
total_samples = 0;

for i = 1:numTest
    % 随机小尺度信道
    h_k = (randn + 1j*randn) / sqrt(2);
    h_mk = (randn + 1j*randn) / sqrt(2);

    % 鲁棒功率分配
    [Pd_robust, Pc_robust] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
        alpha_k_test, alpha_mk_test, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

    % 非鲁棒功率分配
    [Pd_nonrobust, Pc_nonrobust] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
        alpha_k_test, alpha_mk_test, h_k, h_mk, gamma0);

    % 测试中断概率（蒙特卡洛）
    numErr = 1000;
    for j = 1:numErr
        % 生成实际信道误差
        e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
        e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);

        % 实际信道
        hk_actual = epsi_k * h_k + e_k;
        hmk_actual = epsi_mk * h_mk + e_mk;
        gk_actual = alpha_k_test * abs(hk_actual)^2;
        gmk_actual = alpha_mk_test * abs(hmk_actual)^2;

        % 实际SINR
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

fprintf('\n=== 中断概率统计 ===\n');
fprintf('鲁棒算法:   %.4f (%.2f%%)\n', outage_count_robust/total_samples, 100*outage_count_robust/total_samples);
fprintf('非鲁棒算法: %.4f (%.2f%%)\n', outage_count_nonrobust/total_samples, 100*outage_count_nonrobust/total_samples);
fprintf('目标p0:     %.4f (%.2f%%)\n', p0, 100*p0);

%% 详细分析单个案例
fprintf('\n=== 单个案例分析 ===\n');
h_k = (randn + 1j*randn) / sqrt(2);
h_mk = (randn + 1j*randn) / sqrt(2);
fprintf('h_k = %.3f + %.3fj, |h_k|^2 = %.3f\n', real(h_k), imag(h_k), abs(h_k)^2);
fprintf('h_mk = %.3f + %.3fj, |h_mk|^2 = %.3f\n', real(h_mk), imag(h_mk), abs(h_mk)^2);

[Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
    alpha_k_test, alpha_mk_test, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

fprintf('\n分配的功率:\n');
fprintf('Pd_opt = %.4f W (%.1f dBm), Pd_max = %.2f W\n', Pd_opt, 10*log10(Pd_opt*1000), Pd_max);
fprintf('Pc_opt = %.4f W (%.1f dBm), Pc_max = %.2f W\n', Pc_opt, 10*log10(Pc_opt*1000), Pc_max);

% 计算Pc0, Pd0参考值
den0 = alpha_mk_test*(1-epsi_mk^2)*(1/p0-1)*epsi_k^2*abs(h_k)^2 ...
       - (1-epsi_k^2)*alpha_mk_test*epsi_mk^2*abs(h_mk)^2;
Pc0 = (1-epsi_k^2)*sig2/den0;
Pd0 = Pc0*gamma0*alpha_mk_test*(1-epsi_mk^2)*(1-p0)/(alpha_k_test*(1-epsi_k^2)*p0);
fprintf('\nPc0 = %.4f W, Pd0 = %.4f W\n', Pc0, Pd0);
fprintf('Pc_max = %.4f W, Pd_max = %.4f W\n', Pc_max, Pd_max);

% 验证约束条件
B_check = Pd_opt*alpha_k_test*(1-epsi_k^2);
C_check = sig2+Pc_opt*epsi_mk^2*alpha_mk_test*abs(h_mk)^2;
D_check = Pc_opt*alpha_mk_test*(1-epsi_mk^2);
log_tmp_check = -log(1-p0) + epsi_k^2*abs(h_k)^2/(1-epsi_k^2);
log_LHS_check = C_check*gamma0/B_check + log(1+D_check/B_check*gamma0);

fprintf('\n约束验证 (应满足 log_LHS >= log_tmp):\n');
fprintf('log_LHS = %.4f\n', log_LHS_check);
fprintf('log_tmp = %.4f (=-log(1-p0) + epsi^2|h|^2/(1-epsi^2))\n', log_tmp_check);
if log_LHS_check >= log_tmp_check
    status_str = '满足';
else
    status_str = '不满足';
end
fprintf('差值 = %.4f (%s)\n', log_LHS_check - log_tmp_check, status_str);

% 理论预测的中断概率上界
fprintf('\n理论分析:\n');
fprintf('马尔可夫不等式给出的中断概率上界: p0 = %.0e\n', p0);

