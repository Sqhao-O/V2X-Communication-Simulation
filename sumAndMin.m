function [ sumVal, minVal ] = sumAndMin( capMat, assignment )
%SUMANDMIN 计算指派结果中各元素的总和与最小值
%
% 说明：
%   给定一个容量矩阵和一个指派向量（来自 Hungarian 算法），
%   本函数提取被指派元素的容量值，并返回它们的总和与最小值。
%   对于 assignment(i)=0（未指派的行），对应的元素不参与计算。
%
% 输入参数：
%   capMat     - NxM 矩阵，capMat(i,j) 表示将第 j 列指派给第 i 行时的容量值
%   assignment - Nx1 列向量，assignment(i) 表示第 i 行被指派的列索引，
%                 若为 0 则表示该行未指派
%
% 输出参数：
%   sumVal - 标量，所有有效指派元素的容量之和
%   minVal - 标量，所有有效指派元素中的最小容量
%
% 示例：
%   capMat = [1 3 2; 4 2 5; 2 6 1];
%   assignment = [2 0 3]';  % 第1行→第2列，第2行未指派，第3行→第3列
%   [s, m] = sumAndMin(capMat, assignment);  % s=1+1=2, m=1
%
% 注意：
%   若所有行均未指派（assignment(i)=0 对所有 i 成立），minVal 由 min()
%   函数返回 +Inf。
%
% 作者：
%   Le Liang, Georgia Institute of Technology, July 29, 2016

    % 获取矩阵行数和列数
    nrows = size(capMat, 1);
    ncols = size(capMat, 2);

    % 初始化向量，用于存储每个行被指派元素的容量值
    % 未指派的行对应元素保持为 0，不参与后续求和/求最小值
    vec = zeros(nrows, 1);

    % 遍历每一行，提取被指派列的容量值
    for ii = 1 : nrows
        % 取列索引（round 处理浮点精度问题）
        ci = round(assignment(ii));
        if ci >= 1 && ci <= ncols
            % 有效指派：将对应矩阵元素复制到 vec 中
            vec(ii) = capMat(ii, ci);
        else
            % 该行未指派：保持为 0（后续求和时自然排除）
            vec(ii) = 0;
        end
    end

    % 计算聚合指标
    sumVal = sum(vec);   % 所有有效配对的容量总和
    minVal = min(vec);  % 所有有效配对中的最小容量

end
