% DIAGNOSE_OUTAGE2  详细分析鲁棒功率分配算法的内部行为
%
% 目的：
%   在 diagnose_outage.m 的基础上，对单个确定的信道实现进行
%   逐步手工计算，展示算法从输入参数到输出功率的完整流程。
%
% 分析步骤：
%   1. 手动计算参考点 (Pc0, Pd0) — 判断 Case 类型
%   2. 调用 calOptPower — 获取算法输出的 (Pc_opt, Pd_opt)
%   3. 验证 Markov 约束不等式是否满足（检查 lhs >= rhs）
%   4. 蒙特卡洛验证 — 在 10,000 次误差样本上统计真实中断概率
%
% 与 diagnose_outage.m 的区别：
%   - diagnose_outage.m：多信道平均，关注统计趋势
%   - diagnose_outage2.m：单信道逐步分析，关注算法内部机制
%
% 被测函数：calOptPower.m
% -----------------------------------------------------------------

clear; clc;

%% =================================================================
%%  系统参数设置
%% =================================================================
sig2 = 10^(-114/10);     % 噪声功率（线性值，W）
gamma0 = 10^(5/10);      % V2V SINR 阈值（5 dB，线性值）
Pd_max = 10^(23/10);     % V2V 最大功率（23 dBm，线性值）
Pc_max = 10^(23/10);     % V2I 最大功率（23 dBm，线性值）
p0 = 1e-9;               % 目标中断概率（Markov 上界设计值，已收紧）

% 信道时间相关系数（v=60km/h, T=1ms, fc=2GHz）
epsi_k = besselj(0, 2 * pi * (1e-3) * (2e9) * (60/3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('p0 = %.0e, gamma0 = %.2f, epsi_k = %.4f\n', p0, gamma0, epsi_k);

%% =================================================================
%%  选取固定的信道实现（可复现分析）
%% =================================================================
alpha_k = 1e-6;     % V2V 直连链路大尺度衰落因子
alpha_mk = 1e-8;    % V2I→V2V 干扰链路大尺度衰落因子

% 固定小尺度信道（复高斯采样值），确保分析结果可复现
h_k = (-0.898 + 0.828j) / sqrt(2);    % 归一化（|h|^2 ≈ 1.5/2 = 0.75）
h_mk = (0.388 - 0.742j) / sqrt(2);

h_k = h_k * sqrt(2);   % 恢复为单位方差的复高斯（|h|^2 期望 = 1）
h_mk = h_mk * sqrt(2);

fprintf('\n=== 信道状态 ===\n');
fprintf('alpha_k = %.2e, alpha_mk = %.2e\n', alpha_k, alpha_mk);
fprintf('|h_k|^2 = %.3f, |h_mk|^2 = %.3f\n', abs(h_k)^2, abs(h_mk)^2);

%% =================================================================
%%  Step 1：手工计算参考点 (Pc0, Pd0)
%% =================================================================
% (Pc0, Pd0) 是论文中定义的功率空间分界点，用于判断 Case 类型：
%   Case I:  Pd_max <= Pd0               → V2V 功率受限
%   Case II: Pd_max > Pd0 且 Pc_max > Pc0 → V2I 功率受限但在可行域内
%   Case III: 其他                         → 不可行区域（需回退策略）
%
% den0 的公式来自论文中 Pd1 >= Pd0 >= Pd2 的区间定义
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
    fprintf('den0 <= 0! 信道条件极差，(Pc0, Pd0) 无正解\n');
    Pc0 = inf;
    Pd0 = 0;
end

%% =================================================================
%%  Step 2：调用鲁棒功率分配算法
%% =================================================================
[Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

fprintf('\n=== 算法输出 ===\n');
fprintf('Pd_opt = %.4f W (%.1f%% of Pd_max)\n', Pd_opt, 100*Pd_opt/Pd_max);
fprintf('Pc_opt = %.4f W (%.1f%% of Pc_max)\n', Pc_opt, 100*Pc_opt/Pc_max);

%% =================================================================
%%  Step 3：验证 Markov 约束不等式的满足情况
%% =================================================================
% 中间参数定义（与 calOptPower 内部一致）：
%   B = Pd * alpha_k * (1 - epsi_k^2)         — V2V 信号随机分量强度
%   C = sig2 + Pc * alpha_mk * epsi_mk^2 * |h_mk|^2  — 干扰确定分量
%   D = Pc * alpha_mk * (1 - epsi_mk^2)       — 干扰随机分量强度
%
% 约束条件（取对数形式）：
%   C*gamma0/B + log(1 + D*gamma0/B) >= -log(1-p0) + epsi_k^2*|h_k|^2/(1-epsi_k^2)
B = Pd_opt*alpha_k*(1-epsi_k^2);
C = sig2+Pc_opt*epsi_mk^2*alpha_mk*abs(h_mk)^2;
D = Pc_opt*alpha_mk*(1-epsi_mk^2);

log_tmp = -log(1-p0) + epsi_k^2*abs(h_k)^2/(1-epsi_k^2);  % RHS
log_LHS = C*gamma0/B + log(1+D/B*gamma0);                   % LHS

fprintf('\n=== 约束验证 ===\n');
fprintf('  B = %.6e (V2V信号随机分量)\n', B);
fprintf('  C = %.6e (干扰+噪声确定分量)\n', C);
fprintf('  D = %.6e (干扰随机分量)\n', D);
fprintf('  log_LHS = %.4f (当前分配的Markov上界对数)\n', log_LHS);
fprintf('  log_tmp = %.4f (目标约束对数)\n', log_tmp);
fprintf('  差值 = %.4f (%s)\n', log_LHS-log_tmp, iif(log_LHS >= log_tmp, 'OK', 'FAIL'));

%% =================================================================
%%  Step 4：蒙特卡洛验证实际中断概率
%% =================================================================
% 对同一信道估计 (h_k, h_mk)，生成 10,000 次独立的信道误差样本，
% 统计实际 SINR 低于 gamma0 的频率 → 实际中断概率
numMC = 10000;
outage = 0;
for i = 1:numMC
    % 生成信道误差（复高斯，零均值，方差 = 1-epsi^2）
    e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
    e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);

    % 实际信道 = 估计值×相关系数 + 误差项
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
fprintf('目标p0:       %.4f (%.2f%%)\n', p0, 100*p0);
fprintf('实际/目标 比值: %.1e\n', (outage/numMC)/p0);
fprintf('若比值 << 1，说明 Markov 上界远比实际宽松\n');

%% =================================================================
%%  辅助函数：三元条件运算符
%% =================================================================
function s = iif(cond, a, b)
    if cond, s = a; else, s = b; end
end
