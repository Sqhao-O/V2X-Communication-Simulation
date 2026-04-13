% DIAGNOSE_OUTAGE2  详细分析功率分配算法行为

clear; clc;

%% 参数设置
sig2 = 10^(-114/10);
gamma0 = 10^(5/10);
Pd_max = 10^(23/10);
Pc_max = 10^(23/10);
p0 = 1e-9;  % 进一步收紧约束以达到0.1%目标
epsi_k = besselj(0, 2 * pi * (1e-3) * (2e9) * (60/3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('p0 = %.0e, gamma0 = %.2f, epsi_k = %.4f\n', p0, gamma0, epsi_k);

%% 分析单个典型案例
alpha_k = 1e-6;
alpha_mk = 1e-8;
h_k = (-0.898 + 0.828j) / sqrt(2);  % 归一化
h_mk = (0.388 - 0.742j) / sqrt(2);

h_k = h_k * sqrt(2);  % 恢复原始单位
h_mk = h_mk * sqrt(2);

fprintf('\n=== 信道状态 ===\n');
fprintf('alpha_k = %.2e, alpha_mk = %.2e\n', alpha_k, alpha_mk);
fprintf('|h_k|^2 = %.3f, |h_mk|^2 = %.3f\n', abs(h_k)^2, abs(h_mk)^2);

%% 手动计算 Pc0, Pd0
den0 = alpha_mk*(1-epsi_mk^2)*(1/p0-1)*epsi_k^2*abs(h_k)^2 ...
       - (1-epsi_k^2)*alpha_mk*epsi_mk^2*abs(h_mk)^2;
fprintf('\nden0 = %.6e\n', den0);

if den0 > 0
    Pc0 = (1-epsi_k^2)*sig2/den0;
    Pd0 = Pc0*gamma0*alpha_mk*(1-epsi_mk^2)*(1-p0)/(alpha_k*(1-epsi_k^2)*p0);
    fprintf('Pc0 = %.6f W, Pd0 = %.6f W\n', Pc0, Pd0);
    fprintf('Pc_max = %.2f W, Pd_max = %.2f W\n', Pc_max, Pd_max);
    fprintf('Pd_max <= Pd0 ? %s\n', iif(Pd_max <= Pd0, '是 (Case I)', '否'));
    fprintf('Pc_max > Pc0 ? %s\n', iif(Pc_max > Pc0, '是 (Case II)', '否'));
else
    fprintf('den0 <= 0! 信道条件极差\n');
    Pc0 = inf;
    Pd0 = 0;
end

%% 调用算法
[Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

fprintf('\n=== 算法输出 ===\n');
fprintf('Pd_opt = %.4f W (%.1f%% of max)\n', Pd_opt, 100*Pd_opt/Pd_max);
fprintf('Pc_opt = %.4f W (%.1f%% of max)\n', Pc_opt, 100*Pc_opt/Pc_max);

%% 验证约束
B = Pd_opt*alpha_k*(1-epsi_k^2);
C = sig2+Pc_opt*epsi_mk^2*alpha_mk*abs(h_mk)^2;
D = Pc_opt*alpha_mk*(1-epsi_mk^2);
log_tmp = -log(1-p0) + epsi_k^2*abs(h_k)^2/(1-epsi_k^2);
log_LHS = C*gamma0/B + log(1+D/B*gamma0);

fprintf('\n=== 约束验证 ===\n');
fprintf('B = %.6e\n', B);
fprintf('C = %.6e\n', C);
fprintf('D = %.6e\n', D);
fprintf('log_LHS = %.4f\n', log_LHS);
fprintf('log_tmp = %.4f\n', log_tmp);
fprintf('差值 = %.4f (%s)\n', log_LHS-log_tmp, iif(log_LHS >= log_tmp, 'OK', 'FAIL'));

%% 蒙特卡洛验证中断概率
numMC = 10000;
outage = 0;
for i = 1:numMC
    e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
    e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);
    hk_actual = epsi_k * h_k + e_k;
    hmk_actual = epsi_mk * h_mk + e_mk;
    gk_actual = alpha_k * abs(hk_actual)^2;
    gmk_actual = alpha_mk * abs(hmk_actual)^2;
    sinr_actual = Pd_opt * gk_actual / (sig2 + Pc_opt * gmk_actual);
    if sinr_actual < gamma0
        outage = outage + 1;
    end
end

fprintf('\n=== 蒙特卡洛验证 ===\n');
fprintf('样本数: %d\n', numMC);
fprintf('中断次数: %d\n', outage);
fprintf('实际中断概率: %.4f (%.2f%%)\n', outage/numMC, 100*outage/numMC);
fprintf('目标中断概率: %.4f (%.2f%%)\n', p0, 100*p0);

%% 辅助函数
function s = iif(cond, a, b)
    if cond, s = a; else, s = b; end
end
