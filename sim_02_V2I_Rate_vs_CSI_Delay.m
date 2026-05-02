% MAIN_RATE_VS_T  脚本B：V2I 总吞吐量 vs CSI 反馈周期 T
%
% 说明：
%   对比鲁棒算法与非鲁棒算法在不同车速和 CSI 反馈周期下
%   的 V2I 总吞吐量（bps/Hz），验证鲁棒算法对延迟的鲁棒性。
%
%   实验设计：
%     车速列表：v = 50, 100, 150 km/h（3 条曲线）
%     周期列表：T = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 ms（6 个数据点）
%     组合：3 车速 × 6 周期 × 2 算法 = 36 个数据点
%     每点取 channNum 次信道实现的平均值
%
%   预期现象：
%     - 鲁棒算法的吞吐量随 T 增大下降较缓
%     - 非鲁棒算法在 T 大（epsi 小）时吞吐量急剧下降
%
% 输出图表：
%   V2I_Rate_vs_T.png / V2I_Rate_vs_T.pdf
%     - 实线 + 填充标记：鲁棒算法
%     - 虚线 + 空心标记：非鲁棒算法
%     - 蓝色/黑色/绿色：v=50/100/150 km/h
%
% 仿真模式：
%   fastMode = true：channNum = 200（约 17 秒）
%   fastMode = false：channNum = 500（约 3~5 分钟）
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
% fastMode = true：快速验证（~30秒，channNum=200）
% fastMode = false：正式论文结果（~3-5分钟，channNum=500）
fastMode = true;
if fastMode
    channNum = 200;
    fprintf('*** 快速测试模式: channNum=%d ***\n', channNum);
else
    channNum = 1000;
    fprintf('*** 正式仿真模式: channNum=%d ***\n', channNum);
end
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

r0 = 0.5;           % V2I 最低速率要求（bps/Hz）
dB_gamma0 = 5;     % V2V 最低 SINR 阈值（dB）
p0 = 1e-4;         % V2V 目标中断概率（Markov 上界设计值）

%% =================================================================
%%  线性参数转换
%% =================================================================
sig2 = 10^(dB_sig2 / 10);
gamma0 = 10^(dB_gamma0 / 10);
Pd_max = 10^(dB_Pd_max / 10);
Pc_max = 10^(dB_Pc_max / 10);

numCUE = 20;
numDUE = 20;

%% =================================================================
%%  仿真变量配置
%% =================================================================
T_list = [0.2:0.3:0.8, 1.0:0.4:4.6];  % CSI 反馈周期列表（ms），12个采样点更密集
v_list = [50, 100, 150];          % 车速列表（km/h）

% 预分配结果矩阵：v_list × T_list
sumRate_robust = zeros(length(v_list), length(T_list));
sumRate_nonrobust = zeros(length(v_list), length(T_list));

%% =================================================================
%%  主循环：遍历车速和周期组合
%% =================================================================
for v_idx = 1 : length(v_list)
    v = v_list(v_idx);
    d_avg = 2.5 * v / 3.6;  % 平均车头间距（m）
    fprintf('处理车速 v=%d km/h ...\n', v);

    for T_idx = 1 : length(T_list)
        T = T_list(T_idx);

        % 计算时间相关系数（车速和周期共同决定）
        % epsi 越大 → 信道相关性越强 → 延迟误差影响越小
        epsi_k = besselj(0, 2 * pi * (T * 1e-3) * (fc * 1e9) * (v / 3.6) / (3e8));
        epsi_mk = epsi_k;

        % 为当前(v, T)组合设置独立随机种子（避免不同T值间相互影响）
        rng(v_idx * 100 + T_idx);

        % 累加器：对该 (v, T) 组合，多次信道实现求平均
        sumR_robust = 0;
        sumR_nonrobust = 0;
        cntChannLs = 0;

        % -----------------------------------------------------------------
        % Monte Carlo 循环：channNum 次信道实现
        % -----------------------------------------------------------------
        while cntChannLs < channNum

            %% ------ 车辆拓扑生成 ------
            d0 = sqrt(radius^2 - disBstoHwy^2);
            [genFlag, vehPos, indCUE, indDUE, indDUE2] = ...
                genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE);
            if genFlag == 1
                continue;  % 拓扑生成失败，跳过这次
            end

            %% ------ 大尺度衰落计算 ------
            alpha_mB_ = zeros(1, numCUE);
            alpha_k_ = zeros(1, numDUE);
            alpha_kB_ = zeros(1, numDUE);
            alpha_mk_ = zeros(numCUE, numDUE);

            % CUE → 基站 和 CUE → DUE接收端
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

            % DUE发射端 → DUE接收端 和 DUE发射端 → 基站
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

            %% ------ 鲁棒算法 ------
            % 遍历所有 CUE-DUE 对，计算功率分配和 V2I 容量
            C_mk = zeros(numCUE, numDUE);
            Pd_mk_robust = zeros(numCUE, numDUE);
            Pc_mk_robust = zeros(numCUE, numDUE);
            for m = 1 : numCUE
                g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;  % CUE→基站信道增益
                for k = 1 : numDUE
                    g_kB = alpha_kB_(k) * abs(h_kB_(m, k))^2;  % DUE→基站干扰增益

                    % 鲁棒功率分配
                    [Pd_opt, Pc_opt] = calOptPower(1e-6, sig2, Pc_max, Pd_max, ...
                        alpha_k_(k), alpha_mk_(m, k), epsi_k, epsi_mk, ...
                        h_k_(m, k), h_mk_(m, k), p0, gamma0);

                    Pd_mk_robust(m, k) = Pd_opt;
                    Pc_mk_robust(m, k) = Pc_opt;

                    % V2I 容量（bps/Hz）
                    C_mk(m, k) = log2(1 + Pc_opt * g_mB / (sig2 + Pd_opt * g_kB));

                    % V2I 最低速率约束
                    if C_mk(m, k) < r0
                        C_mk(m, k) = -infty;  % 不可行配对
                    end
                end
            end

            % Hungarian 最优配对
            [assign_robust, cost_r] = munkres(-C_mk);
            % 计算有效V2I容量（仅计入V2V未中断的配对）
            valid_r = find(assign_robust > 0);
            sumR_robust_ch = 0;
            for pp = 1 : length(valid_r)
                m_r = valid_r(pp);
                k_r = assign_robust(m_r);
                Pd_r = Pd_mk_robust(m_r, k_r);
                Pc_r = Pc_mk_robust(m_r, k_r);
                % 采样实际V2V信道，检查是否中断
                ek_r = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
                emk_r = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);
                hk_act = epsi_k * h_k_(m_r, k_r) + ek_r;
                hmk_act = epsi_mk * h_mk_(m_r, k_r) + emk_r;
                sinr_act = Pd_r * alpha_k_(k_r) * abs(hk_act)^2 ...
                    / (sig2 + Pc_r * alpha_mk_(m_r, k_r) * abs(hmk_act)^2);
                if sinr_act >= gamma0  % V2V未中断，计入V2I容量
                    g_mB = alpha_mB_(m_r) * abs(h_mB_(m_r))^2;
                    g_kB = alpha_kB_(k_r) * abs(h_kB_(m_r, k_r))^2;
                    cap = log2(1 + Pc_r * g_mB / (sig2 + Pd_r * g_kB));
                    sumR_robust_ch = sumR_robust_ch + cap;
                end
            end
            sumR_robust = sumR_robust + max(0, sumR_robust_ch);

            %% ------ 非鲁棒算法 ------
            C_mk = zeros(numCUE, numDUE);
            Pd_mk_nonrobust = zeros(numCUE, numDUE);
            Pc_mk_nonrobust = zeros(numCUE, numDUE);
            for m = 1 : numCUE
                g_mB = alpha_mB_(m) * abs(h_mB_(m))^2;
                for k = 1 : numDUE
                    g_kB = alpha_kB_(k) * abs(h_kB_(m, k))^2;

                    % 非鲁棒功率分配（假设完美 CSI）
                    [Pd_opt_nr, Pc_opt_nr] = calOptPower_nonrobust(sig2, Pc_max, Pd_max, ...
                        alpha_k_(k), alpha_mk_(m, k), h_k_(m, k), h_mk_(m, k), gamma0);

                    Pd_mk_nonrobust(m, k) = Pd_opt_nr;
                    Pc_mk_nonrobust(m, k) = Pc_opt_nr;

                    C_mk(m, k) = log2(1 + Pc_opt_nr * g_mB / (sig2 + Pd_opt_nr * g_kB));
                    if C_mk(m, k) < r0
                        C_mk(m, k) = -infty;
                    end
                end
            end

            [assign_nonrobust, cost_nr] = munkres(-C_mk);
            % 计算有效V2I容量（仅计入V2V未中断的配对）
            valid_nr = find(assign_nonrobust > 0);
            sumR_nonrobust_ch = 0;
            for pp = 1 : length(valid_nr)
                m_nr = valid_nr(pp);
                k_nr = assign_nonrobust(m_nr);
                Pd_nr = Pd_mk_nonrobust(m_nr, k_nr);
                Pc_nr = Pc_mk_nonrobust(m_nr, k_nr);
                % 采样实际V2V信道，检查是否中断
                ek_nr = sqrt(1 - epsi_k^2) * (randn + 1j*randn) / sqrt(2);
                emk_nr = sqrt(1 - epsi_mk^2) * (randn + 1j*randn) / sqrt(2);
                hk_act_nr = epsi_k * h_k_(m_nr, k_nr) + ek_nr;
                hmk_act_nr = epsi_mk * h_mk_(m_nr, k_nr) + emk_nr;
                sinr_act_nr = Pd_nr * alpha_k_(k_nr) * abs(hk_act_nr)^2 ...
                    / (sig2 + Pc_nr * alpha_mk_(m_nr, k_nr) * abs(hmk_act_nr)^2);
                if sinr_act_nr >= gamma0  % V2V未中断，计入V2I容量
                    g_mB = alpha_mB_(m_nr) * abs(h_mB_(m_nr))^2;
                    g_kB = alpha_kB_(k_nr) * abs(h_kB_(m_nr, k_nr))^2;
                    cap = log2(1 + Pc_nr * g_mB / (sig2 + Pd_nr * g_kB));
                    sumR_nonrobust_ch = sumR_nonrobust_ch + cap;
                end
            end
            sumR_nonrobust = sumR_nonrobust + max(0, sumR_nonrobust_ch);

            cntChannLs = cntChannLs + 1;
        end

        % 信道实现平均：总吞吐量（bps/Hz）
        sumRate_robust(v_idx, T_idx) = sumR_robust / channNum;
        sumRate_nonrobust(v_idx, T_idx) = sumR_nonrobust / channNum;
    end
end

%% =================================================================
%%  保存仿真数据（供 thesis_figures.m 重绘）
%% =================================================================
save(fullfile(outputFolder, 'sim_02_data.mat'), ...
    'T_list', 'v_list', 'sumRate_robust', 'sumRate_nonrobust', 'channNum');

%% =================================================================
%%  绘图（学术规范，单图展示鲁棒vs非鲁棒）
%% =================================================================
LineWidth = 1.5;
MarkerSize = 7;
FontSize = 10.5;
FontName = 'SimSun';  % 宋体（毕设规范）

% 配色方案（学术规范）
colors = [0.0 0.45 0.75; 0.05 0.05 0.05; 0.0 0.6 0.3];  % 蓝、黑、绿
markers = {'o', 's', '^'};  % 圆形、方形、三角形

% 创建单张图（鲁棒vs非鲁棒在同一图中对比）
fig = figure('Position', [100, 80, 1000, 700]);
set(gcf, 'Color', 'white', 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 14], 'PaperPosition', [0.5 0.5 19 13]);

ax = axes('Parent', fig);
set(ax, 'FontName', FontName, 'FontSize', FontSize + 1, 'TickDir', 'in', ...
    'Color', 'white', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5], ...
    'MinorGridAlpha', 0.1, 'MinorGridColor', [0.7 0.7 0.7], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.015 0.015]);
grid(ax, 'on'); grid(ax, 'minor'); hold(ax, 'on');

% 绘制鲁棒算法（实线+填充标记）
h_legend = [];
legend_str = {};
for v_idx = 1 : length(v_list)
    c = colors(v_idx, :);
    mkr = markers{v_idx};
    h = plot(ax, T_list, sumRate_robust(v_idx, :), [mkr, '-'], ...
        'LineWidth', LineWidth, 'MarkerSize', MarkerSize + 1, ...
        'MarkerFaceColor', c, 'Color', c);
    h_legend = [h_legend, h];
    legend_str{end+1} = sprintf('鲁棒, v=%d km/h', v_list(v_idx));
end

% 绘制非鲁棒算法（虚线+空心标记）
for v_idx = 1 : length(v_list)
    c = colors(v_idx, :);
    mkr = markers{v_idx};
    h = plot(ax, T_list, sumRate_nonrobust(v_idx, :), [mkr, '--'], ...
        'LineWidth', LineWidth, 'MarkerSize', MarkerSize + 1, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', c, 'Color', c);
    h_legend = [h_legend, h];
    legend_str{end+1} = sprintf('非鲁棒, v=%d km/h', v_list(v_idx));
end

xlabel('CSI反馈周期 T (ms)', 'FontName', FontName, 'FontSize', FontSize);
ylabel('有效V2I总吞吐量 (bps/Hz)', 'FontName', FontName, 'FontSize', FontSize);
xlim([0, 5]);
ylim([0, max(max(sumRate_robust(:)), max(sumRate_nonrobust(:))) * 1.2]);

legend(ax, h_legend, legend_str, ...
    'FontName', FontName, 'FontSize', 9, ...
    'Location', 'northeast', 'Box', 'off');

hold(ax, 'off');

% 输出图像
setThesisFont(gcf);
print('-dpng', '-r600', fullfile(outputFolder, 'sim_02_V2I_Rate_vs_CSI_Delay.png'));
print('-dpdf', '-r600', fullfile(outputFolder, 'sim_02_V2I_Rate_vs_CSI_Delay.pdf'));
fprintf('图已保存: %s/sim_02_V2I_Rate_vs_CSI_Delay.png / .pdf\n', outputFolder);

%% =================================================================
%%  打印数值结果
%% =================================================================
fprintf('\n===== 有效V2I吞吐量结果 =====\n');
fprintf('说明：有效V2I吞吐量仅计入V2V未中断的配对\n');
disp('非鲁棒算法:');
disp(sumRate_nonrobust);
disp('鲁棒算法:');
disp(sumRate_robust);

toc
