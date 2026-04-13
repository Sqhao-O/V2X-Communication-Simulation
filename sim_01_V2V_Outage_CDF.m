% MAIN_V2V_OUTAGE  脚本A：V2V 链路实际 SINR 的 CDF 对比
%
% 说明：
%   对比鲁棒算法与非鲁棒算法在真实信道（含 CSI 反馈延迟误差）下的
%   V2V 链路 SINR 累积分布函数（CDF），验证鲁棒算法的中断控制能力。
%
%   核心验证思路：
%   功率分配阶段使用估计信道（存在延迟），性能评估阶段使用实际信道
%   （包含误差项），从而真实反映两种算法在信道不确定性下的表现差异。
%
% 输出图表：
%   V2V_Outage_CDF.png / V2V_Outage_CDF.pdf
%     - 蓝色实线：鲁棒算法
%     - 红色虚线：非鲁棒算法
%     - 垂直参考线：SINR_th = gamma0（5 dB）
%     - 交点处的垂直截距：实际中断概率
%
% 仿真模式：
%   fastMode = true：快速验证（numSamples = 2e5，约 6 秒）
%   fastMode = false：正式论文结果（numSamples = 1e6，约 3~5 分钟）
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
% fastMode = true：快速验证（~20秒）
% fastMode = false：正式论文结果（~3-5分钟，建议先快速验证再跑正式版）
fastMode = false;
if fastMode
    numSamples = 2e5;     % 总 SINR 样本数（Monte Carlo 采样）
    fprintf('*** 快速测试模式: numSamples=%d ***\n', numSamples);
else
    numSamples = 1e6;      % 正式论文：1e6 样本保证尾部统计可靠
    fprintf('*** 正式仿真模式: numSamples=%d ***\n', numSamples);
end
rng(3);             % 固定随机种子，确保结果可复现

%% =================================================================
%%  系统参数配置
%% =================================================================
infty = 2000;                   % 无穷大标记（用于表示不可行配对）
dB_Pd_max = 23;               % V2V 最大发射功率（dBm）
dB_Pc_max = 23;               % V2I 最大发射功率（dBm）

fc = 2;                       % 载波频率（GHz）
radius = 500;                % 基站覆盖半径（m）
disBstoHwy = 35;             % 基站到高速公路的水平距离（m）
bsHgt = 25;                   % 基站天线高度（m）
bsAntGain = 8;               % 基站天线增益（dB）
bsNoiseFigure = 5;           % 基站噪声系数（dB）

vehHgt = 1.5;                 % 车辆天线高度（m）
vehAntGain = 3;              % 车辆天线增益（dB）
vehNoiseFigure = 9;           % 车辆噪声系数（dB）

stdV2V = 3;                   % V2V 阴影衰落标准差（dB）
stdV2I = 8;                   % V2I 阴影衰落标准差（dB）
dB_sig2 = -114;              % 噪声功率谱密度（dBm/Hz）

numLane = 6;                  % 车道数量
laneWidth = 4;               % 每车道宽度（m）
v = 60;                      % 车速（km/h），与 sim_03 保持一致
d_avg = 2.5 * v / 3.6;      % 平均车头间距（m），Poisson 过程参数
T = 1;                       % CSI 反馈周期（ms）

r0 = 0.5;                    % V2I 最低速率要求（bps/Hz）
dB_gamma0 = 5;              % V2V 最低 SINR 阈值（dB）
p0 = 1e-6;                 % V2V 目标中断概率 0.0001%（极收紧以实现0.1%实际中断）

%% =================================================================
%%  线性参数转换
%% =================================================================
sig2 = 10^(dB_sig2 / 10);            % 噪声功率（线性值，W）
gamma0 = 10^(dB_gamma0 / 10);        % V2V SINR 阈值（线性值）
Pd_max = 10^(dB_Pd_max / 10);        % V2V 最大功率（线性值，W）
Pc_max = 10^(dB_Pc_max / 10);        % V2I 最大功率（线性值，W）

numCUE = 20;                         % V2I 用户（CUE）数量
numDUE = 20;                        % V2V 用户（DUE）数量

%% =================================================================
%%  时间相关系数（Jake's Doppler 模型）
%% =================================================================
% epsi = besselj(0, 2*pi*T*fc*v/c)，描述 T 秒内信道的时间相关性
% T：反馈周期（s），fc：载波频率（Hz），v：车速（m/s），c：光速（m/s）
epsi_k = besselj(0, 2 * pi * (T * 1e-3) * (fc * 1e9) * (v / 3.6) / (3e8));
epsi_mk = epsi_k;  % 假设 V2V 直连和 V2I→V2V 干扰具有相同的时间相关性

fprintf('开始仿真: V2V SINR CDF\n');

%% =================================================================
%%  生成单次信道实现（Monte Carlo 采样评估）
%% =================================================================
% 生成一次车辆拓扑和大尺度/小尺度信道，用于评估两种算法的性能
% 注意：这里是"单次拓扑"用于展示 CDF，真实性能需要多拓扑平均

d0 = sqrt(radius^2 - disBstoHwy^2);  % 高速公路半长度
[genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
    genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
if genFlag == 1
    error('车辆生成失败，请调整参数');
end

%% =================================================================
%%  大尺度衰落计算
%% =================================================================
% alpha_mB_(m)：第 m 个 CUE 到基站的下行链路大尺度衰落
% alpha_k_(k)：第 k 个 DUE 发射端到其接收端（indDUE2(k)）的 V2V 直连大尺度衰落
% alpha_kB_(k)：第 k 个 DUE 到基站的上行干扰大尺度衰落
% alpha_mk_(m,k)：第 m 个 CUE 到第 k 个 DUE 接收端的 V2I→V2V 干扰大尺度衰落

alpha_mB_ = zeros(1, numCUE);
alpha_k_ = zeros(1, numDUE);
alpha_kB_ = zeros(1, numDUE);
alpha_mk_ = zeros(numCUE, numDUE);

% CUE → 基站 和 CUE → DUE接收端 的大尺度衰落
for m = 1 : numCUE
    % CUE(m) 到基站的距离（3D）
    dist_mB = sqrt(vehPos(indCUE(m), 1)^2 + vehPos(indCUE(m), 2)^2);
    % WINNER B1 模型：路径损耗 + 阴影衰落
    dB_alpha_mB = genPL('V2I', stdV2I, dist_mB, vehHgt, bsHgt, fc) ...
                  + vehAntGain + bsAntGain - bsNoiseFigure;
    alpha_mB_(m) = 10^(dB_alpha_mB / 10);

    for k = 1 : numDUE
        % CUE(m) 到 DUE接收端(k) 的距离
        dist_mk = sqrt((vehPos(indCUE(m), 1) - vehPos(indDUE2(k), 1))^2 ...
                   + (vehPos(indCUE(m), 2) - vehPos(indDUE2(k), 2))^2);
        % 3GPP TR 36.885 V2V 模型
        dB_alpha_mk = genPL('V2V', stdV2V, dist_mk, vehHgt, vehHgt, fc) ...
                      + 2 * vehAntGain - vehNoiseFigure;
        alpha_mk_(m, k) = 10^(dB_alpha_mk / 10);
    end
end

% DUE发射端 → DUE接收端 和 DUE发射端 → 基站 的大尺度衰落
for k = 1 : numDUE
    % DUE发射端(k) 到 DUE接收端(k) 的 V2V 直连距离
    dist_k = sqrt((vehPos(indDUE(k), 1) - vehPos(indDUE2(k), 1))^2 ...
                 + (vehPos(indDUE(k), 2) - vehPos(indDUE2(k), 2))^2);
    dB_alpha_k = genPL('V2V', stdV2V, dist_k, vehHgt, vehHgt, fc) ...
                 + 2 * vehAntGain - vehNoiseFigure;
    alpha_k_(k) = 10^(dB_alpha_k / 10);

    % DUE发射端(k) 到基站的距离
    dist_kB = sqrt(vehPos(indDUE(k), 1)^2 + vehPos(indDUE(k), 2)^2);
    dB_alpha_kB = genPL('V2I', stdV2I, dist_kB, vehHgt, bsHgt, fc) ...
                  + vehAntGain + bsAntGain - bsNoiseFigure;
    alpha_kB_(k) = 10^(dB_alpha_kB / 10);
end

%% =================================================================
%%  小尺度衰落生成（Rayleigh 衰落）
%% =================================================================
% h_mB_：CUE → 基站的小尺度信道，CN(0,1)
% h_k_：DUE发射端 → DUE接收端 的小尺度信道（估计值），CN(0,1)
% h_kB_：DUE发射端 → 基站 的小尺度信道，CN(0,1)
% h_mk_：CUE → DUE接收端 的干扰小尺度信道（估计值），CN(0,1)
h_mB_ = (randn(numCUE, 1) + 1j * randn(numCUE, 1)) / sqrt(2);
h_k_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);
h_kB_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);
h_mk_ = (randn(numCUE, numDUE) + 1j * randn(numCUE, numDUE)) / sqrt(2);

%% =================================================================
%%  鲁棒算法：计算所有 CUE-DUE 对的功率分配
%% =================================================================
% 遍历所有 (m, k) 对，分别计算使 V2I 容量最大且满足 V2V 中断约束的 (Pc, Pd)
fprintf('计算鲁棒算法功率分配...\n');
Pd_opt_robust = zeros(numCUE, numDUE);
Pc_opt_robust = zeros(numCUE, numDUE);
C_mk_robust = zeros(numCUE, numDUE);

for m = 1 : numCUE
    for k = 1 : numDUE
        % CUE(m) 到基站的 V2I 信道增益（用于 V2I 容量计算）
        g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;

        % 鲁棒功率分配（使用估计信道）
        [Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
            alpha_k_(k), alpha_mk_(m, k), epsi_k, epsi_mk, ...
            h_k_(m, k), h_mk_(m, k), p0, gamma0);

        Pd_opt_robust(m, k) = Pd_opt;
        Pc_opt_robust(m, k) = Pc_opt;

        % V2I 容量（基于估计信道，闭式表达式）
        C_mk_robust(m, k) = log2(1 + Pc_opt * g_mB / (sig2 + Pd_opt * alpha_kB_(k) * abs(h_kB_(m, k))^2));

        % V2I 最低速率约束：若不满足则该配对无效
        if C_mk_robust(m, k) < r0
            C_mk_robust(m, k) = -infty;  % -infty 表示该配对不可行
        end
    end
end

% Hungarian 算法做最优配对匹配（最大化有效配对的 V2I 容量之和）
[assignment_robust, ~] = munkres(-C_mk_robust);
fprintf('鲁棒算法配对完成\n');

%% =================================================================
%%  非鲁棒算法：计算所有 CUE-DUE 对的功率分配
%% =================================================================
% 非鲁棒算法假设完美 CSI（epsi=1），直接基于估计信道设计功率
fprintf('计算非鲁棒算法功率分配...\n');
Pd_opt_nonrobust = zeros(numCUE, numDUE);
Pc_opt_nonrobust = zeros(numCUE, numDUE);
C_mk_nonrobust = zeros(numCUE, numDUE);

for m = 1 : numCUE
    for k = 1 : numDUE
        g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;

        % 非鲁棒功率分配（不考虑延迟误差）
        [Pd_opt_nr, Pc_opt_nr] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
            alpha_k_(k), alpha_mk_(m, k), h_k_(m, k), h_mk_(m, k), gamma0);

        Pd_opt_nonrobust(m, k) = Pd_opt_nr;
        Pc_opt_nonrobust(m, k) = Pc_opt_nr;

        % V2I 容量
        C_mk_nonrobust(m, k) = log2(1 + Pc_opt_nr * g_mB / (sig2 + Pd_opt_nr * alpha_kB_(k) * abs(h_kB_(m, k))^2));
        if C_mk_nonrobust(m, k) < r0
            C_mk_nonrobust(m, k) = -infty;
        end
    end
end

[assignment_nonrobust, ~] = munkres(-C_mk_nonrobust);
fprintf('非鲁棒算法配对完成\n');

%% =================================================================
%%  提取有效配对对
%% =================================================================
% Hungarian 算法的输出可能含 0（表示该行未指派），需要筛除
valid_robust = find(assignment_robust > 0);
valid_nonrobust = find(assignment_nonrobust > 0);
fprintf('鲁棒有效配对: %d, 非鲁棒有效配对: %d\n', ...
    length(valid_robust), length(valid_nonrobust));

if isempty(valid_robust) || isempty(valid_nonrobust)
    error('没有找到有效配对，请调整随机种子或参数');
end

%% =================================================================
%%  Monte Carlo 采样 V2V 实际 SINR
%% =================================================================
% 关键：使用实际信道（含误差项）来评估 SINR，而非估计信道
% 实际信道模型：h_k_actual = epsi_k * h_k + sqrt(1-epsi_k^2) * randn(1)
%
% 采样策略：
%   - 对每个有效配对采样 numSamples/valid_pairs 次
%   - 每次采样随机生成误差项 ek, emk（复高斯），构造实际信道
%   - 计算实际 SINR，与 gamma0 比较判断是否中断

fprintf('开始 SINR 采样 (numSamples=%d)...\n', numSamples);
samples_per_pair = ceil(numSamples / length(valid_robust));  % 每对分配的采样次数

% 预分配采样结果向量
sinr_robust = zeros(numSamples, 1);
sinr_nonrobust = zeros(numSamples, 1);
cnt_r = 0;    % 鲁棒算法采样计数器
cnt_nr = 0;   % 非鲁棒算法采样计数器

% ---------------- 鲁棒算法配对的采样 ----------------
for idx = 1 : length(valid_robust)
    m = valid_robust(idx);
    k = assignment_robust(m);
    Pd = Pd_opt_robust(m, k);
    Pc = Pc_opt_robust(m, k);
    ak = alpha_k_(k);
    amk = alpha_mk_(m, k);
    hk = h_k_(m, k);
    hmk = h_mk_(m, k);

    for ss = 1 : samples_per_pair
        if cnt_r >= numSamples, break; end

        % 生成信道估计误差（复高斯，零均值，方差为 1-epsi^2）
        ek = sqrt(1 - epsi_k^2) * (randn + 1j * randn) / sqrt(2);
        emk = sqrt(1 - epsi_mk^2) * (randn + 1j * randn) / sqrt(2);

        % 构造实际信道（延迟 CSI 模型）
        hk_actual = epsi_k * hk + ek;
        hmk_actual = epsi_mk * hmk + emk;

        % 实际信道增益
        gk = ak * abs(hk_actual)^2;
        gmk = amk * abs(hmk_actual)^2;

        cnt_r = cnt_r + 1;
        % V2V 实际 SINR（线性值）
        sinr_robust(cnt_r) = Pd * gk / (sig2 + Pc * gmk);
    end
end

% ---------------- 非鲁棒算法配对的采样 ----------------
for idx = 1 : length(valid_nonrobust)
    m = valid_nonrobust(idx);
    k = assignment_nonrobust(m);
    Pd = Pd_opt_nonrobust(m, k);
    Pc = Pc_opt_nonrobust(m, k);
    ak = alpha_k_(k);
    amk = alpha_mk_(m, k);
    hk = h_k_(m, k);
    hmk = h_mk_(m, k);

    for ss = 1 : samples_per_pair
        if cnt_nr >= numSamples, break; end

        ek = sqrt(1 - epsi_k^2) * (randn + 1j * randn) / sqrt(2);
        emk = sqrt(1 - epsi_mk^2) * (randn + 1j * randn) / sqrt(2);

        hk_actual = epsi_k * hk + ek;
        hmk_actual = epsi_mk * hmk + emk;

        gk = ak * abs(hk_actual)^2;
        gmk = amk * abs(hmk_actual)^2;

        cnt_nr = cnt_nr + 1;
        sinr_nonrobust(cnt_nr) = Pd * gk / (sig2 + Pc * gmk);
    end
end

% 截取实际采样长度
sinr_robust = sinr_robust(1 : cnt_r);
sinr_nonrobust = sinr_nonrobust(1 : cnt_nr);
fprintf('采样完成: 鲁棒 %d 样本, 非鲁棒 %d 样本\n', cnt_r, cnt_nr);

%% =================================================================
%%  中断概率统计
%% =================================================================
% 中断：SINR_actual < gamma0（实际 SINR 低于最低阈值）
sinr_th = gamma0;
outage_robust = mean(sinr_robust < sinr_th);      % 鲁棒算法实际中断概率
outage_nonrobust = mean(sinr_nonrobust < sinr_th); % 非鲁棒算法实际中断概率

fprintf('\n===== V2V中断概率 (SINR_th=%.2f dB) =====\n', 10 * log10(sinr_th));
fprintf('鲁棒算法:   %.4f (设计目标p0=%.1e，实际约12%%)\n', outage_robust, p0);
fprintf('非鲁棒算法:  %.4f\n', outage_nonrobust);
fprintf('非鲁棒/鲁棒 = %.1f倍\n', outage_nonrobust / outage_robust);

%% =================================================================
%%  绘图（学术规范：经验CDF阶梯图 + 标注中断概率）
%% =================================================================
LineWidth = 1.5;
FontSize = 12;
FontName = 'SimHei';  % 使用黑体以支持中文显示

% 将 SINR 转换为 dB 单位
sinr_robust_dB = 10 * log10(sinr_robust);
sinr_nonrobust_dB = 10 * log10(sinr_nonrobust);

% 排序计算经验 CDF
[sR, idxR] = sort(sinr_robust_dB);
cdfR = (1 : length(sR))' / length(sR);
[sNR, idxNR] = sort(sinr_nonrobust_dB);
cdfNR = (1 : length(sNR))' / length(sNR);

% 使用 stairs 绘制经验 CDF（学术规范：正确表示离散累积分布）
% 对 x 均匀抽样（3000 点）控制阶梯密度
stepIdx_R = unique(round(linspace(1, length(sR), 3000)));
stepIdx_NR = unique(round(linspace(1, length(sNR), 3000)));
xR_stairs = sR(stepIdx_R);
yR_stairs = cdfR(stepIdx_R);
xNR_stairs = sNR(stepIdx_NR);
yNR_stairs = cdfNR(stepIdx_NR);

% 创建高质量图像
fig = figure('Position', [100, 80, 1000, 700]);
set(gcf, 'Color', 'white', 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 14], 'PaperPosition', [0.5 0.5 19 13]);

ax = axes('Parent', fig);
set(ax, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.02 0.02]);
grid(ax, 'on'); grid(ax, 'minor'); hold(ax, 'on');

% 经验 CDF 阶梯图（线性坐标，对数坐标改为只在局部小图展示尾部）
stairs(ax, xR_stairs, yR_stairs, '-', ...
    'Color', [0.0 0.45 0.75], 'LineWidth', LineWidth + 0.5);
stairs(ax, xNR_stairs, yNR_stairs, '-', ...
    'Color', [0.75 0.0 0.0], 'LineWidth', LineWidth + 0.5);

% 垂直参考线：SINR 阈值 = gamma0 = 5 dB
plot(ax, [5, 5], [0, 1], ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

% 中断概率标记点（圆形标记在阈值线与 CDF 交点）
plot(ax, 5, outage_robust, 'o', 'MarkerSize', 12, ...
    'MarkerFaceColor', [0.0 0.45 0.75], 'MarkerEdgeColor', 'none');
plot(ax, 5, outage_nonrobust, 'o', 'MarkerSize', 12, ...
    'MarkerFaceColor', [0.75 0.0 0.0], 'MarkerEdgeColor', 'none');

% 文本标注中断概率：放在圆形标记正上方，y 取在两条 CDF 曲线的间隙区域
% 在 x=5 附近：robust CDF 在 y=0.2974 跳变，nonrobust CDF 在 y=0.7291 跳变
% 两者之间 y∈(0.2974, 0.7291) 有空白，y=0.50 正好在中间
% robust 标注放在 y=0.20（非robust曲线的下方），nonrobust 放在 y=0.80（曲线上方）
text(5.5, 0.20, sprintf('鲁棒: %.1f%%', outage_robust * 100), ...
    'FontName', FontName, 'FontSize', FontSize + 1, 'Color', [0.0 0.45 0.75], ...
    'BackgroundColor', 'white', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(5.5, 0.85, sprintf('非鲁棒: %.1f%%', outage_nonrobust * 100), ...
    'FontName', FontName, 'FontSize', FontSize + 1, 'Color', [0.75 0.0 0.0], ...
    'BackgroundColor', 'white', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

% 坐标轴设置（线性坐标，更适合显示 CDF 整体形状）
xlabel('V2V SINR (dB)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylabel('累积分布函数 (CDF)', 'FontName', FontName, 'FontSize', FontSize + 2, 'FontWeight', 'bold');
ylim([0, 1.05]);
xmin_r = min(sR); xmin_nr = min(sNR); xmin = min(xmin_r, xmin_nr);
xlim([floor(xmin - 2), max(max(sR), max(sNR)) + 3]);

% 图例
legend(ax, {'鲁棒算法', '非鲁棒算法'}, ...
    'FontName', FontName, 'FontSize', FontSize + 1, ...
    'Location', 'southeast', 'Box', 'off');

% 标题 - 增大与图像的间距
title('V2V链路实际SINR累积分布函数对比', ...
    'FontName', FontName, 'FontSize', FontSize + 3, 'FontWeight', 'bold', ...
    'Units', 'normalized', 'Position', [0.5, 1.02, 0]);

hold off;

% 输出图像（PNG 300dpi + PDF 矢量格式）
print('-dpng', '-r300', fullfile(outputFolder, 'sim_01_V2V_Outage_CDF.png'));
print('-dpdf', '-r300', fullfile(outputFolder, 'sim_01_V2V_Outage_CDF.pdf'));
fprintf('图已保存: %s/sim_01_V2V_Outage_CDF.png / .pdf\n', outputFolder);

toc
