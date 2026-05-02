function setThesisFont(figHandle)
% SETTHESISFONT  按本科毕设规范设置图中所有文字的字体
%   中文 -> SimSun (宋体)，英文/数字 -> Times New Roman
%   坐标轴标签/刻度/标题：10.5pt (五号字)
%   图例/标注文字：9pt (小五号)

    if nargin < 1 || ~isvalid(figHandle)
        figHandle = gcf;
    end

    % ---- 1. 全局默认字体 ----
    set(figHandle, 'DefaultAxesFontName', 'Times New Roman', ...
                   'DefaultAxesFontSize', 10.5);

    % ---- 2. 所有坐标轴 ----
    allAxes = findall(figHandle, 'Type', 'axes');
    labelHandles = [];  % 收集标签句柄，供第4步排除
    for k = 1 : length(allAxes)
        ax = allAxes(k);
        % 刻度字体：Times New Roman 10.5pt + Box封闭 + 刻度朝内
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 10.5, ...
            'Box', 'on', 'TickDir', 'in', 'LineWidth', 0.5);

        % XLabel
        xlh = get(ax, 'XLabel');
        if isvalid(xlh)
            setLabelFont(xlh);
            labelHandles = [labelHandles; xlh];
        end

        % YLabel
        ylh = get(ax, 'YLabel');
        if isvalid(ylh)
            setLabelFont(ylh);
            labelHandles = [labelHandles; ylh];
        end

        % Title
        th = get(ax, 'Title');
        if isvalid(th)
            setLabelFont(th);
            labelHandles = [labelHandles; th];
        end
    end

    % ---- 3. 图例：统一 SimSun 9pt（图例含中英文混合）----
    allLegends = findall(figHandle, 'Type', 'legend');
    for k = 1 : length(allLegends)
        set(allLegends(k), 'FontName', 'SimSun', 'FontSize', 9);
    end

    % ---- 4. 所有 text 对象（排除坐标轴标签和标题）----
    allTexts = findall(figHandle, 'Type', 'text');
    for k = 1 : length(allTexts)
        obj = allTexts(k);
        if ~isempty(labelHandles) && any(labelHandles == obj)
            continue;  % 标签/标题已在第2步处理，跳过
        end
        str = get(obj, 'String');
        if iscell(str)
            str = strjoin(str, ' ');
        end
        if hasChinese(str)
            set(obj, 'FontName', 'SimSun', 'FontSize', 9);
        else
            set(obj, 'FontName', 'Times New Roman', 'FontSize', 9);
        end
    end

    % ---- 5. sgtitle (R2018b+) ----
    try
        sgTitles = findobj(figHandle, 'Type', 'subplottext');
        for k = 1 : length(sgTitles)
            set(sgTitles(k), 'FontName', 'SimSun', 'FontSize', 10.5);
        end
    catch ME
        % 旧版本无 sgtitle 对象，忽略
    end

    % ---- 6. colorbar ----
    allCbars = findall(figHandle, 'Type', 'colorbar');
    for k = 1 : length(allCbars)
        set(allCbars(k), 'FontName', 'Times New Roman', 'FontSize', 10.5);
    end
end

%% 辅助函数：根据字符串内容设置 label 字体
function setLabelFont(h)
    str = get(h, 'String');
    if iscell(str)
        str = strjoin(str, ' ');
    end
    if hasChinese(str)
        set(h, 'FontName', 'SimSun', 'FontSize', 10.5);
    else
        set(h, 'FontName', 'Times New Roman', 'FontSize', 10.5);
    end
end

%% 辅助函数：检测字符串是否包含中文
function flag = hasChinese(str)
    flag = false;
    if isempty(str) || ~ischar(str) && ~isstring(str)
        return;
    end
    str = char(str);
    % 检查是否有 CJK 统一表意文字 (U+4E00-U+9FFF)
    flag = ~isempty(regexp(str, '[\u4E00-\u9FFF]', 'once'));
end
