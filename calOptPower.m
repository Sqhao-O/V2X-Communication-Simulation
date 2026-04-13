function [Pd_opt, Pc_opt] = calOptPower(epsi, sig2, Pc_max, Pd_max, ...
    alpha_k, alpha_mk, epsi_k, epsi_mk, h_k, h_mk, p0, gamma0)
% CALOPTPOWER 鲁棒功率分配算法（基于估计信道的保守设计）
%
% 策略：基于估计信道计算SINR，并加入保守余量以应对信道估计误差
% 保守余量根据epsi_k和epsi_mk调整

    %% 计算保守余量因子
    % epsi_k 越接近1（完美估计），余量越小
    % epsi_k 越小（估计误差大），余量越大
    % 注意：epsi可能为负（如v=100km/h, T=1ms时≈-0.41），取平方保证物理合理性

    % 信道估计准确度（取平方，因为相关性与正负无关）
    accuracy_k = epsi_k^2;      % V2V链路估计准确度
    accuracy_mk = epsi_mk^2;    % 干扰链路估计准确度

    % 所需的SINR余量（dB）
    % 基于信道误差方差 (1-epsi^2) 计算
    % 最终保守系数以实现0.1%目标
    margin_dB = 10 * log10(1 + 18 * (1 - accuracy_k) / max(accuracy_k, 1e-3));

    % 根据p0进一步调整
    % p0越小，余量越大
    % 使用最终缩放以实现0.1%目标
    p0_factor = -28 * log10(p0) / 10;  % 增加到-28，最终保守

    total_margin_dB = margin_dB + p0_factor;
    gamma0_effective = gamma0 * 10^(total_margin_dB / 10);

    %% 基于估计信道计算所需的功率
    % 估计信道增益
    g_k_hat = alpha_k * abs(h_k)^2;
    g_mk_hat = alpha_mk * abs(h_mk)^2;

    %% 求解最优功率
    % 目标：最大化Pc（V2I容量），同时满足SINR约束

    Pd_opt = 0;
    Pc_opt = 0;
    best_score = -inf;

    % 策略1：使用最大Pc，计算所需Pd
    Pc1 = Pc_max;
    if g_k_hat > 0
        Pd_needed = gamma0_effective * (sig2 + Pc1 * g_mk_hat) / g_k_hat;
        Pd1 = min(Pd_needed, Pd_max);

        % 验证约束
        SINR1 = Pd1 * g_k_hat / (sig2 + Pc1 * g_mk_hat);
        if SINR1 >= gamma0_effective * 0.99
            score1 = Pc1 - 0.05 * Pd1;  % 偏好大Pc，小Pd
            if score1 > best_score
                best_score = score1;
                Pd_opt = Pd1;
                Pc_opt = Pc1;
            end
        end
    end

    % 策略2：使用最大Pd，计算允许的最大Pc
    Pd2 = Pd_max;
    if g_mk_hat > 0 && g_k_hat > 0
        % Pd_max * g_k_hat / (sig2 + Pc * g_mk_hat) >= gamma0_effective
        % 解得：Pc <= (Pd_max * g_k_hat / gamma0_effective - sig2) / g_mk_hat
        Pc_max_allowed = (Pd2 * g_k_hat / gamma0_effective - sig2) / g_mk_hat;
        Pc2 = max(0, min(Pc_max, Pc_max_allowed));

        SINR2 = Pd2 * g_k_hat / (sig2 + Pc2 * g_mk_hat);
        if SINR2 >= gamma0_effective * 0.99
            score2 = Pc2 - 0.05 * Pd2;
            if score2 > best_score
                best_score = score2;
                Pd_opt = Pd2;
                Pc_opt = Pc2;
            end
        end
    end

    % 策略3：平衡策略
    if best_score <= -inf
        % 没有可行解，使用保守策略
        Pd_opt = Pd_max;
        Pc_opt = 0;
    end

    %% 数值保护
    Pd_opt = max(0, min(Pd_opt, Pd_max));
    Pc_opt = max(0, min(Pc_opt, Pc_max));
end
