function [assignment,cost] = munkres(costMat)
% MUNKRES   Munkres (Hungarian) Algorithm for Linear Assignment Problem.
%
% 在论文中的作用：
%   用于 V2I 和 V2V 用户配对问题的最优匹配。给定 numCUE × numDUE
%   的 V2I 容量矩阵 C_mk（或 -C_mk 转化为最小化问题），Hungarian
%   算法找出使总 V2I 容量最大化的 CUE-DUE 一一配对方案。
%
%   不可行配对（V2I 容量 < r0 或 V2V 中断约束不满足）的代价设
%   为 -infty（最大化问题），转换后为 +infty（最小化问题）。
%   Hungarian 算法会避开这些不可行配对，确保最终匹配中的每对
%   都具有可行的功率分配。
%
%   调用方式（论文中的典型用法）：
%     [assignment, ~] = munkres(-C_mk);
%     % assignment(m) = k 表示 CUE m 配对给 DUE k
%     % assignment(m) = 0 表示 CUE m 未配对（无可行 DUE）
%
% [ASSIGN,COST] = MUNKRES(COSTMAT) returns the optimal column indices,
% ASSIGN assigned to each row and the minimum COST based on the assignment
% problem represented by the COSTMAT, where the (i,j)th element represents
% the cost to assign the jth column to the ith row.
%
% Partial assignment: supports rectangular matrices and partial assignments
% where some rows have no valid column (encoded as Inf). For a partial
% assignment, zero elements in the returning ASSIGN vector indicate
% unassigned rows. The cost only contains the cost of assigned tasks.
%
% EXAMPLES:
%   [assignment,cost] = munkres(magic(5)); % 5x5 assignment
%   A=rand(10,7); A(A>0.7)=Inf; [a,b]=munkres(A); % partial assignment
%
% Reference:
%   "Munkres' Assignment Algorithm, Modified for Rectangular Matrices",
%   http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
%
% Version 2.3 by Yi Cao at Cranfield University on 11th September 2011

% ---- 初始化 ----
assignment = zeros(1,size(costMat,1));
cost = 0;

% 标记有效元素（非 NaN 且非 Inf）
validMat = costMat == costMat & costMat < Inf;
% 用一个极大值替代无效元素，确保算法不会选中它们
bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
costMat(~validMat) = bigM;

% 提取有效行列
validCol = any(validMat,1);
validRow = any(validMat,2);

nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return  % 无有效元素，直接返回
end

% 构造方阵（不足部分用极大值填充）
maxv=10*max(costMat(validMat));
dMat = zeros(n) + maxv;
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);

% ---- 代价矩阵的行列标准化（加速收敛）----
% 等价于经典的 "行减最小值、列减最小值" 两步
minR = min(dMat,[],2);                          % 每行的最小值
minC = min(bsxfun(@minus, dMat, minR));         % 行减后每列的最小值

% 零元素的位置矩阵 zP（标记 dMat(i,j) == minR(i) + minC(j) 的位置）
zP = dMat == bsxfun(@plus, minC, minR);

% 初始匹配：贪心选零（无冲突的行列对）
starZ = zeros(n,1);   % starZ(r)=c 表示行 r 的星号零在第 c 列
while any(zP(:))
    [r,c]=find(zP,1);   % 取第一个零
    starZ(r)=c;          % 标记为星号零
    zP(r,:)=false;       % 该行不能再选
    zP(:,c)=false;       % 该列不能再选
end

while 1
    % STEP 3: 检查覆盖 — 若每列都有星号零，则匹配已最大化
    if all(starZ>0)
        break
    end
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;             % 有星号零的列被覆盖
    coverRow = false(n,1);                         % 行覆盖初始化为无
    primeZ = zeros(n,1);                           % primeZ(r)=c 表示行 r 的画线零
    % 找未覆盖区域中的零
    [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
    while 1
        % STEP 4: 在未覆盖区域找零，画线
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);           % 未覆盖零所在行
            uZc = cIdx(1);           % 未覆盖零所在列
            primeZ(uZr) = uZc;       % 标记为画线零
            stz = starZ(uZr);        % 该行是否有星号零
            if ~stz
                Step = 5;            % 该行无星号零 → 进入 STEP 5 增广
                break;
            end
            coverRow(uZr) = true;    % 覆盖该行
            coverColumn(stz) = false;% 取消覆盖星号零所在列
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            % 在新暴露的列中找零
            z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            % STEP 6: 调整 minR/minC 以创造新零
            [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));
            minC(~coverColumn) = minC(~coverColumn) + minval;  % 未覆盖列加 minval
            minR(coverRow) = minR(coverRow) - minval;          % 覆盖行减 minval
        else
            break
        end
    end
    % STEP 5: 交替路径增广 — 沿星号零和画线零交替翻转
    rowZ1 = find(starZ==uZc);       % 找 uZc 列的星号零所在行
    starZ(uZr)=uZc;                 % 将画线零升级为星号零
    while ~isempty(rowZ1)
        starZ(rowZ1)=0;             % 取消原有星号零
        uZc = primeZ(rowZ1);        % 沿画线零前进
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);   % 找下一个星号零
        starZ(uZr)=uZc;             % 将画线零升级为星号零
    end
end

% ---- 将结果映射回原始矩阵的行列索引 ----
rowIdx = find(validRow);     % 有效行在原矩阵中的索引
colIdx = find(validCol);     % 有效列在原矩阵中的索引
starZ = starZ(1:nRows);      % 截取有效部分的星号零
vIdx = starZ <= nCols;       % 星号零对应有效列的行
assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));  % 填入原索引
% 清理无效指派（对应原矩阵中为 Inf 的配对）
pass = assignment(assignment>0);
pass(~diag(validMat(assignment>0,pass))) = 0;
assignment(assignment>0) = pass;
% 计算总代价（trace = 对角线之和 = 各指派对的代价总和）
cost = trace(costMat(assignment>0,assignment(assignment>0)));

% ---- 辅助函数：计算未覆盖区域的最小值 ----
function [minval,rIdx,cIdx]=outerplus(M,x,y)
ny=size(M,2);
minval=inf;
for c=1:ny
    M(:,c)=M(:,c)-(x+y(c));
    minval = min(minval,min(M(:,c)));
end
[rIdx,cIdx]=find(M==minval);
