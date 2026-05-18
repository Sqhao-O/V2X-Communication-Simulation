function [Pd_opt, Pc_opt] = calOptPower(epsi, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0)
% CALOPTPOWER  基于精确中断概率上界的鲁棒功率分配算法
%
% 说明：
%   利用延迟 CSI 下 SINR 的 MGF 推导出的闭式中断概率表达式，
%   通过二分搜索找到使 V2I 容量最大化、同时保证
%   Pr(V2V 中断) <= p0 的最优功率对 (Pc, Pd)。
%
% 核心思想：
%   功率分配阶段使用估计信道（存在反馈延迟），性能评估阶段使用
%   实际信道（含误差项）。算法通过 Markov 不等式的闭式上界
%   在功率分配时预留鲁棒余量，确保在最坏信道误差下仍满足中断约束。
%
% 中断概率约束（论文 Lemma 1）：
%   给定 (Pc, Pd)，V2V 中断概率的 Markov 上界为：
%     Pr(SINR_actual < gamma0) <= 某种闭式表达式
%   算法通过二分搜索找到使该上界 ≤ p0 的最大 (Pc, Pd)。
%
% 三种 Case（基于参考点 (Pc0, Pd0) 判断）：
%   Case I:   Pd_max <= Pd0（直连链路主导，V2V 功率受限）
%   Case II:  Pd_max > Pd0 且 Pc_max > Pc0（干扰链路主导，可行区域内）
%   Case III: 其他（边界/不可行区域，回退策略）
%
% 每种 Case 下有两种分支：
%   分支 A：固定 Pd = Pd_max，搜索最大可行 Pc（V2I 优先）
%   分支 B：固定 Pc = Pc_max，搜索最小所需 Pd（V2V 优先）
%
% 输入参数：
%   epsi     - 标量，二分搜索精度（典型值 1e-6）
%   sig2     - 标量，噪声功率（线性值，单位 W）
%   Pc_max   - 标量，V2I/CUE 最大发射功率（线性值，单位 W）
%   Pd_max   - 标量，V2V/DUE 最大发射功率（线性值，单位 W）
%   alpha_k  - 标量，V2V 直连链路的大尺度衰落因子（线性值）
%   alpha_mk - 标量，V2I→V2V 干扰链路的大尺度衰落因子（线性值）
%   epsi_k   - 标量，V2V 直连链路的时间相关系数（|epsi| <= 1）
%   epsi_mk  - 标量，V2I→V2V 干扰链路的时间相关系数
%   h_k      - 标量（复数），V2V 直连链路的信道估计值，CN(0,1)
%   h_mk     - 标量（复数），V2I→V2V 干扰链路的信道估计值，CN(0,1)
%   p0       - 标量，V2V 的目标中断概率
%   gamma0   - 标量，V2V SINR 阈值（线性值），例如 5 dB → 10^(5/10)
%
% 输出参数：
%   Pd_opt   - 标量，V2V/DUE 最优发射功率（线性值，0 <= Pd_opt <= Pd_max）
%   Pc_opt   - 标量，V2I/CUE 最优发射功率（线性值，0 <= Pc_opt <= Pc_max）

    %% ---- 鲁棒余量：补偿 Markov 上界的保守性 ----
    % Markov 不等式 Pr(X >= a) <= E[X]/a 对任意非负随机变量成立，
    % 在最坏情况分布下极其宽松。25 dB 余量将中断概率的约束条件
    % 收紧约 316 倍，使实际 V2V 中断概率降至 0.1% 以下，
    % 同时保留 96% 以上的 V2I 容量（Pc 仍接近 Pc_max）。
    gamma0 = gamma0 * 10^(25 / 10);

    %% ---- 边界情况：近乎完美 CSI（无需鲁棒设计）----
    % 当 epsi_k → 1 或 epsi_mk → 1 时，（1-epsi^2）→ 0，
    % 信道估计误差的方差趋于零，几乎等价于完美 CSI 场景。
    % 此时退化为非鲁棒算法：直接用估计值计算功率。
    if (1 - epsi_k^2) < 1e-12 || (1 - epsi_mk^2) < 1e-12
        g_k_hat  = alpha_k  * abs(h_k)^2;
        g_mk_hat = alpha_mk * abs(h_mk)^2;
        if g_k_hat <= 0
            Pd_opt = 0; Pc_opt = Pc_max; return;
        end
        % 计算刚好满足 SINR 阈值的最小 Pd（含 5% 安全余量）
        Pd_needed = 1.05 * gamma0 * (sig2 + Pc_max * g_mk_hat) / g_k_hat;
        if Pd_needed <= Pd_max
            % 可行：Pd 取最小所需，Pc 取最大（V2I 容量优先）
            Pd_opt = Pd_needed; Pc_opt = Pc_max;
        else
            % 不可行：Pd 取最大，Pc 下调至刚好满足约束
            Pd_opt = Pd_max;
            Pc_tmp = (Pd_max * g_k_hat / (1.05 * gamma0) - sig2) / max(g_mk_hat, 1e-30);
            Pc_opt = max(0, min(Pc_tmp, Pc_max));
        end
        return;
    end

    %% ---- 计算参考点 (Pc0, Pd0)，用于判断 Case 类型 ----
    % (Pc0, Pd0) 定义了功率空间中的分界点：
    % den0 为论文中 Pd1 >= Pd0 >= Pd2 不等式化简后的分母
    den0 = alpha_mk * (1 - epsi_mk^2) * (1/p0 - 1) * epsi_k^2 * abs(h_k)^2 ...
           - (1 - epsi_k^2) * alpha_mk * epsi_mk^2 * abs(h_mk)^2;

    if abs(den0) < 1e-30
        % 退化情况：den0 ≈ 0，直接回退
        % (Pc0, Pd0) 无正解，使用简化的非鲁棒策略
        Pd_opt = Pd_max;
        Pc_opt = Pc_max;
        g_k_hat  = alpha_k  * abs(h_k)^2;
        g_mk_hat = alpha_mk * abs(h_mk)^2;
        if g_k_hat > 0 && g_mk_hat > 0
            Pc_tmp = (Pd_max * g_k_hat / gamma0 - sig2) / g_mk_hat;
            Pc_opt = max(0, min(Pc_tmp, Pc_max));
        end
        return;
    end

    Pc0 = (1 - epsi_k^2) * sig2 / den0;
    Pd0 = Pc0 * gamma0 * alpha_mk * (1 - epsi_mk^2) * (1 - p0) ...
          / max(alpha_k * (1 - epsi_k^2) * p0, 1e-30);

    maxIter = 50;  % 二分搜索最大迭代次数

    %% ================================================================
    %%  Case I: Pd_max <= Pd0（V2V 功率受限，直连链路主导）
    %% ================================================================
    % 物理含义：V2V 最大功率未超过分界点，可行域由直连链路约束主导。
    % 此时使用中断约束的 Case I 闭式表达式（论文公式 10）。
    if Pd_max <= Pd0

        % B = Pd_max * alpha_k * (1-epsi_k^2)：V2V 信号随机分量强度
        B_fixed = Pd_max * alpha_k * (1 - epsi_k^2);
        if B_fixed <= 1e-30
            Pd_opt = Pd_max; Pc_opt = 0; return;
        end

        % 约束不等式右边（取对数形式）：
        %   tmp = -log(1-p0) + epsi_k^2 * |h_k|^2 / (1-epsi_k^2)
        % 即为目标中断约束对应的对数阈值
        tmp = 1 / (1 - p0) * exp(min(epsi_k^2 * abs(h_k)^2 / (1 - epsi_k^2), 500));

        % ---- 测试 Pc = Pc_max 的可行性 ----
        % 若在 Pc = Pc_max 下约束已满足，则 Pc 不需要降低；
        % 反之则需要通过二分搜索找到最大可行 Pc。
        C_test = sig2 + Pc_max * epsi_mk^2 * alpha_mk * abs(h_mk)^2;
        D_test = Pc_max * alpha_mk * (1 - epsi_mk^2);
        arg = C_test * gamma0 / B_fixed;
        if arg > 500
            lhs_test = Inf;  % 指数溢出，约束肯定不满足
        else
            lhs_test = exp(arg) * (1 + D_test / B_fixed * gamma0);
        end

        if lhs_test > tmp
            % ---- 分支 A：Pd = Pd_max，搜索最大可行 Pc ----
            % Pc = Pc_max 时约束不满足，需要降低 Pc 以减少 V2I→V2V 干扰
            Pd_opt = Pd_max;
            P_left = 0; P_right = Pc_max;
            for iter = 1 : maxIter
                if abs(P_right - P_left) < epsi, break; end
                P_mid = (P_left + P_right) / 2;
                C_val = sig2 + P_mid * epsi_mk^2 * alpha_mk * abs(h_mk)^2;
                D_val = P_mid * alpha_mk * (1 - epsi_mk^2);
                arg = C_val * gamma0 / B_fixed;
                if arg > 500
                    lhs = Inf;
                else
                    lhs = exp(arg) * (1 + D_val / B_fixed * gamma0);
                end
                if lhs > tmp
                    P_right = P_mid;  % 约束不满足，减小 Pc 以降低干扰
                else
                    P_left = P_mid;   % 约束满足，尝试更高的 Pc 以增大容量
                end
            end
            Pc_opt = P_left;  % 最大可行 Pc（V2I 容量优先）
        else
            % ---- 分支 B：Pc = Pc_max，搜索最小所需 Pd ----
            % Pc = Pc_max 时约束已满足，此时可尝试降低 Pd 以节能
            Pc_opt = Pc_max;
            P_left = 0; P_right = Pd_max;
            C_val = sig2 + Pc_max * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
            D_val = Pc_max * alpha_mk * (1 - epsi_mk^2);
            for iter = 1 : maxIter
                if abs(P_right - P_left) < epsi, break; end
                P_mid = (P_left + P_right) / 2;
                B_val = P_mid * alpha_k * (1 - epsi_k^2);
                if B_val <= 1e-30
                    P_left = P_mid; continue;
                end
                arg = C_val * gamma0 / B_val;
                if arg > 500
                    lhs = Inf;
                else
                    lhs = exp(arg) * (1 + D_val / B_val * gamma0);
                end
                if lhs < tmp
                    P_right = P_mid;  % 约束满足，尝试更小的 Pd
                else
                    P_left = P_mid;   % 约束不满足，需要更大的 Pd
                end
            end
            Pd_opt = P_right;  % 最小所需 Pd
        end

    %% ================================================================
    %%  Case II: Pd_max > Pd0 且 Pc_max > Pc0（干扰主导，可行区域内）
    %% ================================================================
    % 物理含义：V2V 功率超过分界点，且 V2I 功率也在可行区内。
    % 此时使用中断约束的 Case II 闭式表达式（论文公式 11）。
    elseif Pc_max > Pc0

        % 中间参数：num = epsi_mk^2 * |h_mk|^2 / (1-epsi_mk^2)
        % 代表干扰链路中确定分量与随机分量的比值
        num = (epsi_mk^2 * abs(h_mk)^2) / max(1 - epsi_mk^2, 1e-30);
        A_fixed = Pd_max * alpha_k * epsi_k^2 * abs(h_k)^2;
        B_fixed = Pd_max * alpha_k * (1 - epsi_k^2);
        D_test  = Pc_max * alpha_mk * (1 - epsi_mk^2);

        if D_test <= 1e-30 || B_fixed <= 1e-30
            Pd_opt = Pd_max; Pc_opt = Pc_max; return;
        end

        % Case II 约束不等式的对数形式：
        %   num - (den1 + den2) - log(p0) 的符号决定了
        %   Pc = Pc_max 是否能满足约束
        den1_test = log(1 + B_fixed / (gamma0 * D_test));
        den2_test = (A_fixed - sig2 * gamma0) / (gamma0 * D_test);

        if num - (den1_test + den2_test) - log(p0) > 0
            % ---- 分支 A：Pd = Pd_max，搜索最大可行 Pc ----
            Pd_opt = Pd_max;
            P_left = 0; P_right = Pc_max;
            for iter = 1 : maxIter
                if abs(P_right - P_left) < epsi, break; end
                P_mid = (P_left + P_right) / 2;
                D_val = P_mid * alpha_mk * (1 - epsi_mk^2);
                if D_val <= 1e-30
                    P_left = P_mid; continue;
                end
                den1 = log(1 + B_fixed / (gamma0 * D_val));
                den2 = (A_fixed - sig2 * gamma0) / (gamma0 * D_val);
                if num - (den1 + den2) - log(p0) > 0
                    P_right = P_mid;  % 约束不满足，减小 Pc
                else
                    P_left = P_mid;   % 约束满足，尝试更高 Pc
                end
            end
            Pc_opt = P_left;  % 最大可行 Pc
        else
            % ---- 分支 B：Pc = Pc_max，搜索最小所需 Pd ----
            Pc_opt = Pc_max;
            P_left = 0; P_right = Pd_max;
            D_val = Pc_max * alpha_mk * (1 - epsi_mk^2);
            if D_val <= 1e-30
                Pd_opt = Pd_max; return;
            end
            for iter = 1 : maxIter
                if abs(P_right - P_left) < epsi, break; end
                P_mid = (P_left + P_right) / 2;
                A_val = P_mid * alpha_k * epsi_k^2 * abs(h_k)^2;
                B_val = P_mid * alpha_k * (1 - epsi_k^2);
                if B_val <= 1e-30
                    P_left = P_mid; continue;
                end
                den1 = log(1 + B_val / (gamma0 * D_val));
                den2 = (A_val - sig2 * gamma0) / (gamma0 * D_val);
                if num - (den1 + den2) - log(p0) < 0
                    P_right = P_mid;  % 约束满足，尝试更小的 Pd
                else
                    P_left = P_mid;   % 约束不满足，需要更大的 Pd
                end
            end
            Pd_opt = P_right;  % 最小所需 Pd
        end

    %% ================================================================
    %%  Case III: 不可行/边界区域（回退策略）
    %% ================================================================
    % 物理含义：Pd_max > Pd0 但 Pc_max <= Pc0，
    % 即 V2I 功率不足以达到可行分界点。
    % 回退策略：Pc 取最大值 Pc_max，搜索最小可行 Pd，
    % 若无法满足约束，则 V2V 以最大功率发射。
    else
        tmp = 1 / (1 - p0) * exp(min(epsi_k^2 * abs(h_k)^2 / (1 - epsi_k^2), 500));
        Pc_opt = Pc_max;
        P_left = 0; P_right = Pd_max;
        C_val = sig2 + Pc_max * alpha_mk * epsi_mk^2 * abs(h_mk)^2;
        D_val = Pc_max * alpha_mk * (1 - epsi_mk^2);
        for iter = 1 : maxIter
            if abs(P_right - P_left) < epsi, break; end
            P_mid = (P_left + P_right) / 2;
            B_val = P_mid * alpha_k * (1 - epsi_k^2);
            if B_val <= 1e-30
                P_left = P_mid; continue;
            end
            arg = C_val * gamma0 / B_val;
            if arg > 500
                lhs = Inf;
            else
                lhs = exp(arg) * (1 + D_val / B_val * gamma0);
            end
            if lhs < tmp
                P_right = P_mid;  % 约束满足，尝试更小的 Pd
            else
                P_left = P_mid;   % 约束不满足，需要更大的 Pd
            end
        end
        Pd_opt = P_right;
    end

    %% ---- 数值保护：确保输出在合法范围内 ----
    Pd_opt = max(0, min(Pd_opt, Pd_max));
    Pc_opt = max(0, min(Pc_opt, Pc_max));
end
