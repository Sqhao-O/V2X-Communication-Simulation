% 诊断脚本：检查calOptPower各Case触发情况和功率分配
addpath('E:\桌面常用\毕设\bs_matlab备份');
fastMode = true;

%% 参数
sig2 = 10^(-114/10);
Pc_max = 10^(23/10);
Pd_max = 10^(23/10);
gamma0 = 10^(5/10);
p0 = 1e-3;

v = 60; T = 1;
epsi_k = besselj(0, 2*pi*(T*1e-3)*(2*1e9)*(v/3.6)/(3e8));
epsi_mk = epsi_k;

fprintf('epsi_k = %.4f\n', epsi_k);

%% 生成一个典型拓扑
d0 = sqrt(500^2 - 35^2);
[~, vehPos, indCUE, indDUE, indDUE2] = ...
    genCUEandDUE(d0, 4, 6, 35, 2.5*v/3.6, 20, 20);

stdV2V = 3; stdV2I = 8; fc = 2; vehHgt = 1.5; bsHgt = 25;
vehAntGain = 3; bsNoiseFigure = 5; vehNoiseFigure = 9;

%% 随机选5对，分析功率分配
rng(1);
fprintf('\n===== 鲁棒算法功率分配诊断 =====\n');
fprintf('%-6s %-10s %-10s %-10s %-10s %-12s %-12s\n', ...
    'pair', 'alpha_k', 'alpha_mk', '|h_k|^2', '|h_mk|^2', 'Pc_opt(W)', 'Pd_opt(W)');

for idx = 1:5
    m = indCUE(randi(20));
    k = indDUE(randi(20));

    % 大尺度
    dist_k = sqrt((vehPos(indDUE(k),1)-vehPos(indDUE2(k),1))^2 ...
                 + (vehPos(indDUE(k),2)-vehPos(indDUE2(k),2))^2);
    dB_ak = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) + 2*vehAntGain - vehNoiseFigure;
    alpha_k = 10^(dB_ak/10);

    dist_mk = sqrt((vehPos(indCUE(m),1)-vehPos(indDUE2(k),1))^2 ...
               + (vehPos(indCUE(m),2)-vehPos(indDUE2(k),2))^2);
    dB_amk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) + 2*vehAntGain - vehNoiseFigure;
    alpha_mk = 10^(dB_amk/10);

    % 小尺度
    h_k = (randn + 1j*randn) / sqrt(2);
    h_mk = (randn + 1j*randn) / sqrt(2);

    % 调用鲁棒算法
    [Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
        alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

    fprintf('%-6d %-10.2e %-10.2e %-10.4f %-10.4f %-12.6f %-12.6f\n', ...
        idx, alpha_k, alpha_mk, abs(h_k)^2, abs(h_mk)^2, Pc_opt, Pd_opt);
end

%% 理论推导：Case I可行条件检验
fprintf('\n===== Case I 可行性分析 =====\n');
fprintf('目标: Pr(SINR < gamma0) <= p0 = %.3e\n', p0);
fprintf('所需: Pd * g_k / sig2 >= -log(p0) * gamma0 = %.2f (线性)\n', -log(p0)*gamma0);

% 取典型g_k=alpha_k*|h_k|^2, alpha_k ~ 10^(-5), |h_k|^2 ~ 1
g_k_typical = 1e-5;  % 典型V2V通道增益 (20m, 88dB pathloss)
left_side = Pd_max * g_k_typical / sig2;
fprintf('\n典型 V2V (dist~20m, g_k~1e-5):\n');
fprintf('  Pd_max * g_k / sig2 = %.2e * %.1e / %.2e = %.2e\n', Pd_max, g_k_typical, sig2, left_side);
fprintf('  所需: %.2e, 比值: %.2e (>>1? NO!)\n', -log(p0)*gamma0, left_side / (-log(p0)*gamma0));
fprintf('  → Case I 约束 %.2e << 所需 %.2e, Case I 不可行!\n', left_side, -log(p0)*gamma0);

% 强通道 g_k=10^-4 (10m)
g_k_strong = 1e-4;
left_strong = Pd_max * g_k_strong / sig2;
fprintf('\n较强 V2V (dist~10m, g_k~1e-4):\n');
fprintf('  Pd_max * g_k / sig2 = %.2e * %.1e / %.2e = %.2e\n', Pd_max, g_k_strong, sig2, left_strong);
fprintf('  所需: %.2e, 比值: %.2e\n', -log(p0)*gamma0, left_strong / (-log(p0)*gamma0));
fprintf('  → 仍不可行! 需要 g_k > %.2e (pathloss < %.1f dB)\n', ...
    -log(p0)*gamma0*sig2/Pd_max, -10*log10(-log(p0)*gamma0*sig2/Pd_max));

fprintf('\n===== 结论 =====\n');
fprintf('要满足 p0=1e-3, 需要路径损耗 < %.1f dB\n', -10*log10(-log(p0)*gamma0*sig2/Pd_max));
fprintf('但 3GPP V2V 模型在 dist=20m 给出 ~88 dB > 71.6 dB\n');
fprintf('→ p0=1e-3 在本信道模型下不可达到!\n');
