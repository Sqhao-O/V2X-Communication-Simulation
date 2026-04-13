function [ combinedPL ] = genPL(linkType, stdShadow, dist, hgtTX, hgtRX, freq)
% GENPL  计算大尺度衰落：路径损耗 + 对数正态阴影
%
% 说明：
%   计算 V2I 或 V2V 链路的总大尺度路径损耗（dB 单位），并叠加对数正态
%   阴影衰落。函数返回的是 dB 单位下的总衰落因子（内部已取负号）。
%
%   实现了两种传播模型：
%     - V2I：WINNER B1 模型（郊区宏蜂窝场景）
%     - V2V：3GPP TR 36.885 模型（高速公路场景，近场衍射区）
%
% 输入参数：
%   linkType   - 字符串（不区分大小写）：'V2I' 或 'V2V'
%   stdShadow  - 标量，对数正态阴影的标准差（dB）
%                典型值：stdV2I = 8 dB，stdV2V = 3 dB
%   dist       - 标量，发射天线与接收天线之间的水平距离（米）
%   hgtTX      - 标量，发射天线高度（米）
%   hgtRX      - 标量，接收天线高度（米）
%   freq       - 标量，载波频率（GHz）
%
% 输出参数：
%   combinedPL - 标量，dB 单位下的总大尺度衰落值。
%                内部已取负号，因此路径损耗越大（PL_dB 越正），
%                combinedPL 越负，信道增益越小。
%                转换为线性尺度：alpha = 10^(combinedPL/10)
%
% 路径损耗模型：
%
%   【V2I — WINNER B1，郊区宏蜂窝】
%     PL_dB = 128.1 + 37.6 * log10( sqrt((hgtTX - hgtRX)^2 + dist^2) / 1000 )
%     其中 dist_3D = sqrt((hgtTX - hgtRX)^2 + dist^2) 为三维距离（千米）。
%
%   【V2V — 3GPP TR 36.885，高速公路场景，Table 6.2-1】
%     断点距离（远场-近场转换点）：
%       d_bp = 4 * (hgtTX - 1) * (hgtRX - 1) * freq * 10^9 / (3 * 10^8)
%     近场区（dist <= 3 m）：          PL = 22.7*log10(3) + 41.0 + 20*log10(freq/5)
%     近场衍射区（3 < dist <= d_bp）：PL = 22.7*log10(dist) + 41.0 + 20*log10(freq/5)
%     远场自由空间（dist > d_bp）：    PL = 40*log10(dist) + 9.45
%                                        - 17.3*log10((hgtTX-1)*(hgtRX-1))
%                                        + 2.7*log10(freq/5)
%
% 使用示例：
%   % V2I 链路：基站（h=25m）到车辆（h=1.5m），距离 500m
%   PL = genPL('V2I', 8, 500, 25, 1.5, 2.0);
%   alpha = 10^(PL/10);  % 转换为线性尺度的信道增益因子
%
%   % V2V 链路：两辆车间距 20m，天线高度均为 1.5m
%   PL = genPL('V2V', 3, 20, 1.5, 1.5, 2.0);  % 使用 3GPP 高速公路模型
%
% 作者：
%   Le Liang, Georgia Institute of Technology, 2016年7月10日

    % -----------------------------------------------------------------
    %  V2V 路径损耗模型（3GPP TR 36.885，高速公路场景）
    % -----------------------------------------------------------------
    if strcmp(upper(linkType), 'V2V')
        % 计算断点距离（远场-近场转换点）
        % d_bp = 4 * (h_tx - 1) * (h_rx - 1) * f_c * c^(-1)
        % 其中 c = 3e8 m/s，频率单位为 Hz（故乘以 10^9）
        d_bp = 4*(hgtTX-1)*(hgtRX-1)*freq*10^9/(3*10^8);

        % 模型系数（3GPP TR 36.885，Table 6.2-1）
        A = 22.7;  % 近场区对数斜率（dB）
        B = 41.0;  % 截距（dB）
        C = 20;    % 频率依赖系数（dB）

        % 根据距离分段计算路径损耗
        if dist <= 3
            % 极近距（< 3 m）：使用 3m 处的固定值，
            % 避免 log10(dist) 在 dist → 0 时出现负无穷
            PL = A*log10(3) + B + C*log10(freq/5);
        elseif dist <= d_bp
            % 近场衍射区（3m ~ 断点距离）：对数距离路径损耗模型
            PL = A*log10(dist) + B + C*log10(freq/5);
        else
            % 远场自由空间区（> 断点距离）：频率依赖自由空间 + 地面反射
            PL = 40*log10(dist)+9.45-17.3*log10((hgtTX-1)*(hgtRX-1))+2.7*log10(freq/5);
        end

    % -----------------------------------------------------------------
    %  V2I 路径损耗模型（WINNER B1，郊区宏蜂窝）
    % -----------------------------------------------------------------
    else
        % WINNER B1 模型：适用于郊区宏蜂窝场景（基站对车辆）
        % 适用于天线高度 > 1m、距离可达数公里的场景
        % 采用三维距离：sqrt((hgtTX - hgtRX)^2 + dist^2)，单位转换为千米
        %
        % PL_dB = 128.1 + 37.6 * log10( d_3D_km )
        % 其中 d_3D_km = sqrt((hgtTX - hgtRX)^2 + dist^2) / 1000
        %
        PL = 128.1 + 37.6*log10(sqrt((hgtTX-hgtRX)^2+dist^2)/1000);
    end

    % -----------------------------------------------------------------
    %  叠加对数正态阴影衰落，并取负号
    % -----------------------------------------------------------------
    % 阴影衰落：X_sigma ~ N(0, sigma^2)，建模为 randn * stdShadow
    % 总衰落：-(PL_dB + X_sigma)
    %
    % 取负号的含义：
    %   - 路径损耗越大（PL_dB 越正）→ combinedPL 越负 → alpha 越小
    %   - 正的阴影（最坏情况）→ combinedPL 更负 → alpha 更小
    % 转换为线性尺度的信道增益因子：alpha = 10^(combinedPL / 10)
    % -----------------------------------------------------------------
    combinedPL = - (randn(1)*stdShadow + PL);

end
