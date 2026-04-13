function [Flag, vehPos, indCUE, indDUE, indDUE2] = ...
    genCUEandDUE(d0, laneWidth, numLane, disBstoHwy, d_avg, numCUE, numDUE)
% GENCUEANDDUE  基于空间泊松过程生成高速公路车辆拓扑
%
% 说明：
%   首先按照空间泊松过程（由平均车头间距 d_avg 决定密度）在高速公路
%   每条车道上生成车辆位置；然后从中随机选取 V2I 用户（CUE）和
%   V2V 用户（DUE）发射端，并为每个 DUE 发射端分配其对应的 DUE 接收端
%   （距离最近的另一辆车）。
%
%   本函数生成的是同一时刻的车辆静态拓扑，用于一次信道实现的仿真。
%   多次调用可模拟不同信道实现（蒙特卡洛仿真）。
%
% 输入参数：
%   d0         - 标量，高速公路半长度（米），车辆在
%                [-d0, +d0] 范围内生成（水平方向）
%   laneWidth   - 标量，每条车道的宽度（米）
%   numLane     - 标量，车道数量
%   disBstoHwy  - 标量，基站到高速公路的水平距离（米）
%   d_avg       - 标量，平均车头间距（米），d_avg = 2.5 * v / 3.6，
%                其中 v 为车速（km/h）。泊松过程的强度参数为 1/d_avg。
%   numCUE      - 标量，V2I 用户（CUE）的数量
%   numDUE      - 标量，V2V 用户（DUE）发射端的数量
%
% 输出参数：
%   Flag       - 标量，0 = 成功生成，1 = 生成失败（车辆数不足）
%   vehPos     - Nx2 矩阵，所有车辆的位置坐标（单位：米）
%                第1列：水平坐标（沿高速公路方向，范围 [-d0, d0]）
%                第2列：垂直坐标（垂直于高速公路方向，即到基站的横向距离）
%   indCUE     - numCUEx1 列向量，vehPos 中 CUE 车辆的索引
%   indDUE     - numDUEx1 列向量，vehPos 中 DUE 发射端的索引
%   indDUE2    - numDUEx1 列向量，vehPos 中 DUE 接收端的索引，
%                indDUE2(i) 为 indDUE(i) 所对应的 DUE 接收端索引
%
% 算法流程：
%   Step 1（车辆生成）：
%        每条车道独立按照泊松过程生成车辆，强度为 2*d0 / d_avg
%        （长度为 2*d0 的路段上期望车辆数）。
%        水平坐标在 [-d0, d0] 上均匀分布，纵向坐标固定为
%        disBstoHwy + lane*laneWidth（各车道中心线位置）。
%   Step 2（DUE 选取）：
%        从所有车辆中随机无放回抽取 numDUE 个作为 DUE 发射端。
%   Step 3（DUE 接收端匹配）：
%        对每个 DUE 发射端，在剩余车辆中找到与其欧氏距离最近的车辆
%        作为其对应的 DUE 接收端（保证 V2V 直连链路长度尽可能短）。
%   Step 4（CUE 选取）：
%        从剩余车辆中随机选取 numCUE 个作为 CUE（V2I 用户）。
%
% 示例：
%   % 生成 20 个 CUE 和 20 个 DUE（v=60km/h，T=1ms）
%   d0 = sqrt(500^2 - 35^2);  % 高速公路半长
%   v = 60; d_avg = 2.5 * v / 3.6;
%   [Flag, vehPos, indCUE, indDUE, indDUE2] = ...
%       genCUEandDUE(d0, 4, 6, 35, d_avg, 20, 20);
%   if Flag == 0, disp('生成成功'); end
%
% 注意：
%   - 生成的车辆数必须至少为 numCUE + 2*numDUE，
%     否则函数返回 Flag=1（生成失败）。
%   - 所有坐标单位为米（m）。
%   - 同一辆车不能同时作为 CUE 和 DUE。
%
% 作者：
%   Le Liang, Georgia Institute of Technology, 2017年1月25日

    % -----------------------------------------------------------------
    % 初始化输出变量
    % -----------------------------------------------------------------
    vehPos = [];   % 所有车辆位置
    indCUE = [];   % CUE 索引
    indDUE = [];   % DUE 发射端索引
    indDUE2 = [];  % DUE 接收端索引
    Flag = 0;      % 默认生成成功

    %% =================================================================
    %%  Step 1：生成所有车辆位置（空间泊松过程）
    %% =================================================================
    % 每条车道独立按泊松过程生成车辆
    % 泊松过程强度：lambda = 1/d_avg（辆/米）
    % 在长度为 2*d0 的路段上，期望车辆数为 2*d0 / d_avg
    for ilane = 1 : numLane
        % 本车道的期望车辆数（泊松随机变量）
        npoints = poissrnd(2 * d0 / d_avg);
        % 水平坐标：在 [-d0, d0] 上均匀分布
        pproc = (rand(npoints, 1) * 2 - 1) * d0;
        % 纵向坐标：各车道中心线位置 = disBstoHwy + ilane * laneWidth
        pproc2 = [pproc, (disBstoHwy + ilane * laneWidth) * ones(length(pproc), 1)];
        % 追加到总车辆列表
        vehPos = [vehPos; pproc2];
    end

    % 获取总车辆数
    numVeh = size(vehPos, 1);

    % 检查车辆数是否足够（CUE + DUE发射端 + DUE接收端）
    if numVeh < numCUE + 2 * numDUE
        Flag = 1;  % 生成失败
        return;
    end

    %% =================================================================
    %%  Step 2：随机选取 DUE 发射端
    %% =================================================================
    % 对所有车辆索引做随机排列
    indPerm = randperm(numVeh);
    % 取前 numDUE 个作为 DUE 发射端
    indDUE = indPerm(1 : numDUE);
    % 为每个 DUE 发射端分配 DUE 接收端（见 Step 3）
    indDUE2 = zeros(1, numDUE);

    %% =================================================================
    %%  Step 3：为每个 DUE 发射端分配最近的 DUE 接收端
    %% =================================================================
    % 对于每个 DUE 发射端 ii，在所有车辆中找与之距离最近且尚未被使用的车辆
    % （不能是其他 DUE 发射端或已被分配的 DUE 接收端）
    for ii = 1 : numDUE
        minDist = 2 * d0;  % 初始化为最大可能距离
        tmpInd = 0;        % 暂存找到的最近车辆索引

        % 遍历所有车辆
        for iii = 1 : numVeh
            % 排除自身（iii in indDUE）及已被指派的接收端（iii in indDUE2）
            if any(abs(iii - indDUE) < 1e-6) || any(abs(iii - indDUE2) < 1e-6)
                continue;  % 跳过不满足条件的车辆
            end

            % 计算 DUE(ii) 与车辆 iii 之间的欧氏距离
            newDist = sqrt((vehPos(indDUE(ii), 1) - vehPos(iii, 1))^2 ...
                         + (vehPos(indDUE(ii), 2) - vehPos(iii, 2))^2);

            % 更新最近距离和对应索引
            if newDist < minDist
                tmpInd = iii;
                minDist = newDist;
            end
        end

        % 将找到的最近车辆指派为 DUE(ii) 的接收端
        indDUE2(ii) = tmpInd;
    end

    %% =================================================================
    %%  Step 4：随机选取 CUE（V2I 用户）
    %% =================================================================
    % CUE 必须从剩余车辆中选取（即不在 indDUE 也不在 indDUE2 中）
    cntCUE = numDUE + 1;  % indPerm 中前 numDUE 个已被用作 DUE
    while cntCUE <= numVeh
        if any(abs(indPerm(cntCUE) - indDUE2) < 1e-6)
            % 该车辆已被某个 DUE 接收端使用，跳过
        else
            % 该车辆可用，加入 CUE 列表
            indCUE = [indCUE, indPerm(cntCUE)];
        end
        cntCUE = cntCUE + 1;

        % 达到所需 CUE 数量后退出
        if length(indCUE) >= numCUE
            break;
        end
    end

    % 转换为列向量（确保维度一致）
    indCUE = indCUE(:);
    indDUE = indDUE(:);
    indDUE2 = indDUE2(:);

end
