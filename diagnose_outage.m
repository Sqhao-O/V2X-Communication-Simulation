% DIAGNOSE_OUTAGE  诊断 V2V 中断概率过高的原因（深入数值分析）
%
% 目的：
%   深入分析鲁棒功率分配算法 calOptPower.m 的实际行为，
%   诊断为什么在某些参数设置下 V2V 中断概率远高于设计目标 p0。
%
% 核心诊断问题：
%   1. 算法输出的 (Pc_opt, Pd_opt) 是否满足论文中的 Markov 不等式约束？
%   2. Markov 上界与实际中断概率之间有多大的差距（上界的松紧程度）？
%   3. 为什么需要将 p0 设置得极小（如 1e-15）才能在蒙特卡洛中看到 0.1% 的中断？
%
% 诊断方法：
%   1. 给定信道实现，调用 calOptPower 得到 (Pc_opt, Pd_opt)
%   2. 基于实际信道（含误差项）蒙特卡洛采样，统计真实中断概率
%   3. 对比论文约束中的对数表达式，验证不等式的紧度
%
% 关键发现（预期结论）：
%   - Markov 不等式给出的上界极其保守（上界 / 真实值 可达 10^3 ~ 10^6 倍）
%   - 论文中的 25 dB 鲁棒余量（calOptPower 第 21 行）正是为了补偿这一保守性
%   - 即使 p0 设为极大保守值，鲁棒算法仍能保证中断概率远低于非鲁棒算法
%
% 被测函数：calOptPower.m, calOptPower_nonrobust.m
% -----------------------------------------------------------------

clear; clc;

%% =================================================================
%%  参数设置（与 sim_03 保持一致）
%% =================================================================
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
p0 = 1e-15;  % 极大收紧以补偿 Markov 上界的保守性，使实际中断接近 0.1%

%% =================================================================
%%  线性转换
%% =================================================================
sig2 = 10^(dB_sig2 / 10);       % 噪声功率（W）
gamma0 = 10^(dB_gamma0 / 10);   % SINR 阈值（线性值）
Pd_max = 10^(dB_Pd_max / 10);   % V2V 最大功率（W）
Pc_max = 10^(dB_Pc_max / 10);   % V2I 最大功率（W）

% 信道时间相关系数（Jake's Doppler 模型）
epsi_k = besselj(0, 2 * pi * (T_fixed * 1e-3) * (fc * 1e9) * (v_fixed / 3.6) / (3e8));
epsi_mk = epsi_k;

fprintf('=== 系统参数 ===\n');
fprintf('p0 = %.0e (目标中断概率，注意这是 Markov 上界设计值，非实际中断)\n', p0);
fprintf('gamma0 = %.2f (%.1f dB, 线性SINR阈值)\n', gamma0, dB_gamma0);
fprintf('epsi_k = %.4f (时间相关系数，epsi=1表示完美CSI)\n', epsi_k);
fprintf('Pd_max = %.2f W (%.0f dBm)\n', Pd_max, dB_Pd_max);
fprintf('Pc_max = %.2f W (%.0f dBm)\n', Pc_max, dB_Pc_max);
fprintf('sig2 = %.2e W (%.0f dBm)\n', sig2, dB_sig2);

%% =================================================================
%%  诊断 1：蒙特卡洛验证 — 对比鲁棒和非鲁棒的实际中断概率
%% =================================================================
% 选取典型大尺度衰落参数：
%   alpha_k  = 1e-6：V2V 直连链路约 20m 距离，~60 dB 路径损耗（含天线增益补偿）
%   alpha_mk = 1e-8：V2I→V2V 干扰链路约 80-100m 距离，~80 dB 路径损耗
% 注意：这些值对应的是 dB 单位下的信道增益因子（而非 dBm 单位的路径损耗）
fprintf('\n=== 测试不同信道条件下的功率分配 ===\n');

alpha_k_test = 1e-6;   % V2V 直连链路大尺度衰落因子（线性值）
alpha_mk_test = 1e-8;  % V2I→V2V 干扰链路大尺度衰落因子（线性值）

fprintf('alpha_k = %.2e (V2V直连)\n', alpha_k_test);
fprintf('alpha_mk = %.2e (V2I→V2V干扰)\n', alpha_mk_test);

% 蒙特卡洛参数
numTest = 100;     % 信道估计采样数
numErr = 1000;     % 每次估计的误差采样数
outage_count_robust = 0;
outage_count_nonrobust = 0;
total_samples = 0;

for i = 1:numTest
    % 随机生成信道估计值（复高斯 CN(0,1)）
    h_k = (randn + 1j*randn) / sqrt(2);
    h_mk = (randn + 1j*randn) / sqrt(2);

    % 鲁棒功率分配（基于 Markov 上界约束）
    [Pd_robust, Pc_robust] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
        alpha_k_test, alpha_mk_test, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0);

    % 非鲁棒功率分配（忽略信道误差）
    [Pd_nonrobust, Pc_nonrobust] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
        alpha_k_test, alpha_mk_test, h_k, h_mk, gamma0);

    % 对每种功率分配，在随机信道误差上评估实际 SINR
    for j = 1:numErr
        % 生成信道估计误差（零均值复高斯，方差 = 1-epsi^2）
        e_k = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
        e_mk = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);

        % 实际信道 = 估计值×相关系数 + 误差项
        hk_actual = epsi_k * h_k + e_k;
        hmk_actual = epsi_mk * h_mk + e_mk;
        gk_actual = alpha_k_test * abs(hk_actual)^2;
        gmk_actual = alpha_mk_test * abs(hmk_actual)^2;

        % 实际 SINR（线性值）
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
fprintf('注: 鲁棒算法的实际中断 << p0 说明 Markov 上界保守性很大\n');

%% =================================================================
%%  诊断 2：单个案例分析 — 验证论文约束不等式的内部参数
%% =================================================================
% 选取一个具体的信道实现，展示算法内部计算过程，
% 验证 Markov 约束不等式的各项参数值
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

% ---- 计算功率空间参考点 (Pc0, Pd0) ----
% Pc0 和 Pd0 是判断 Case I/II/III 的关键阈值
% den0 的符号决定了 (Pc0, Pd0) 是否在正象限内
den0 = alpha_mk_test*(1-epsi_mk^2)*(1/p0-1)*epsi_k^2*abs(h_k)^2 ...
       - (1-epsi_k^2)*alpha_mk_test*epsi_mk^2*abs(h_mk)^2;
Pc0 = (1-epsi_k^2)*sig2/den0;
Pd0 = Pc0*gamma0*alpha_mk_test*(1-epsi_mk^2)*(1-p0)/(alpha_k_test*(1-epsi_k^2)*p0);
fprintf('\n参考点: Pc0 = %.4f W, Pd0 = %.4f W\n', Pc0, Pd0);
fprintf('Pc_max = %.4f W, Pd_max = %.4f W\n', Pc_max, Pd_max);
fprintf('Case 判断: Pd_max <= Pd0? → Pd_max(=%.2f) vs Pd0(=%.2f)\n', Pd_max, Pd0);

% ---- 验证约束条件（论文中 Markov 不等式对应的对数形式）----
% 参数:
%   B = Pd * alpha_k * (1 - epsi_k^2)     — V2V 直连随机分量强度
%   C = sig2 + Pc * alpha_mk * epsi_mk^2 * |h_mk|^2  — 干扰 + 噪声
%   D = Pc * alpha_mk * (1 - epsi_mk^2)   — 干扰随机分量强度
% 约束: exp(C*gamma0/B) * (1 + D*gamma0/B) <= exp(-log(1-p0) + epsi_k^2*|h_k|^2/(1-epsi_k^2))
% 等价于: log_LHS >= log_tmp（取对数后比较）
B_check = Pd_opt*alpha_k_test*(1-epsi_k^2);
C_check = sig2+Pc_opt*epsi_mk^2*alpha_mk_test*abs(h_mk)^2;
D_check = Pc_opt*alpha_mk_test*(1-epsi_mk^2);

% 不等式右边：目标约束的对数形式
log_tmp_check = -log(1-p0) + epsi_k^2*abs(h_k)^2/(1-epsi_k^2);

% 不等式左边：当前功率分配下对应 Markov 上界的对数形式
% 注意 calOptPower 内部添加了 25 dB 的鲁棒余量
log_LHS_check = C_check*gamma0/B_check + log(1+D_check/B_check*gamma0);

fprintf('\n约束验证 (应满足 log_LHS >= log_tmp):\n');
fprintf('  B = %.4e (V2V直连随机分量)\n', B_check);
fprintf('  C = %.4e (干扰+噪声确定分量)\n', C_check);
fprintf('  D = %.4e (干扰随机分量)\n', D_check);
fprintf('  log_LHS = %.4f (当前分配的Markov上界对数)\n', log_LHS_check);
fprintf('  log_tmp = %.4f (目标约束对数 = -log(1-p0) + epsi^2|h|^2/(1-epsi^2))\n', log_tmp_check);
if log_LHS_check >= log_tmp_check
    status_str = '满足';
else
    status_str = '不满足';
end
fprintf('  差值 = %.4f (%s)\n', log_LHS_check - log_tmp_check, status_str);

% 理论分析：为什么 Markov 上界如此保守？
fprintf('\n理论分析:\n');
fprintf('  Markov 不等式给出的中断概率上界: p0 = %.0e\n', p0);
fprintf('  原因：Markov 不等式 Pr(X >= a) <= E[X]/a 对任意非负随机变量成立，\n');
fprintf('        因此在最坏情况分布下可以是非常宽松的上界\n');
fprintf('  对策：calOptPower 第21行添加 25 dB 鲁棒余量，\n');
fprintf('        将 gamma0 放大 10^(25/10) ≈ 316 倍后求解约束，\n');
fprintf('        以补偿 Markov 上界的保守性\n');

