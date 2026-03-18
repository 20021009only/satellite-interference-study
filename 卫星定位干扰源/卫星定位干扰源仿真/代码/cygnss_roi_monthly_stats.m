%% cygnss_roi_monthly_stats.m
clear; clc;

%% ===================== 1) 路径与文件 =====================
dataDir = 'D:\CYGNSS_L1_v32\2021\2021.01\';   % 改成你的数据目录
files1 = dir(fullfile(dataDir, '**', '*.nc'));
files2 = dir(fullfile(dataDir, '**', '*.nc4'));
fileList = [files1; files2];

if isempty(fileList)
    error('没有找到 .nc 或 .nc4 文件，请检查 dataDir');
end

%% ===================== 2) ROI 定义 =====================
% 这里把边界替换成你之前已经验证过的 ROI 范围
ROIs = struct( ...
    'name',   {'MiddleEast_NorthAfrica','EastAsia_NWPacific','SouthIndianOcean','SouthPacific'}, ...
    'latMin', {NaN, NaN, NaN, NaN}, ...
    'latMax', {NaN, NaN, NaN, NaN}, ...
    'lonMin', {NaN, NaN, NaN, NaN}, ...
    'lonMax', {NaN, NaN, NaN, NaN}  ...
);

for i = 1:numel(ROIs)
    if any(isnan([ROIs(i).latMin, ROIs(i).latMax, ROIs(i).lonMin, ROIs(i).lonMax]))
        error('请先把 4 个 ROI 的经纬度边界补全');
    end
end

%% ===================== 3) 参数 =====================
minAntCount = 50;   % 计算 Port-Starboard 差值时，每侧至少要有这么多样本
rec = struct( ...
    'roi', {}, ...
    'month', {}, ...
    'noise_sum', {}, 'noise_n', {}, ...
    'rfi_sum', {},   'rfi_n', {}, ...
    'kurt_sum', {},  'kurt_n', {}, ...
    'total_n', {}, ...
    'port_noise_sum', {}, 'port_n', {}, ...
    'star_noise_sum', {}, 'star_n', {} ...
);

recIdx = 0;

%% ===================== 4) 主循环 =====================
for k = 1:numel(fileList)
    fn = fullfile(fileList(k).folder, fileList(k).name);

    fileDate = parseDateFromFilename(fileList(k).name);
    if isnat(fileDate)
        warning('文件名中未识别出日期，跳过: %s', fileList(k).name);
        continue;
    end
    monthDT = dateshift(fileDate, 'start', 'month');

    try
        noise = double(ncread(fn, 'ddm_noise_floor'));
        rfi   = double(ncread(fn, 'rfi_detected'));
        kurt  = double(ncread(fn, 'ddm_kurtosis'));
        ant   = double(ncread(fn, 'ddm_ant'));
        lat   = double(ncread(fn, 'sp_lat'));
        lon   = double(ncread(fn, 'sp_lon'));
    catch ME
        warning('读取失败，跳过文件: %s\n原因: %s', fileList(k).name, ME.message);
        continue;
    end

    % 拉平成列向量
    noise = noise(:);
    rfi   = rfi(:);
    kurt  = kurt(:);
    ant   = ant(:);
    lat   = lat(:);
    lon   = lon(:);

    % 长度一致性检查
    N = min([numel(noise), numel(rfi), numel(kurt), numel(ant), numel(lat), numel(lon)]);
    noise = noise(1:N);
    rfi   = rfi(1:N);
    kurt  = kurt(1:N);
    ant   = ant(1:N);
    lat   = lat(1:N);
    lon   = lon(1:N);

    %% ===== 关键预处理：你这轮已经确认过的正确写法 =====
    noise(noise <= 0) = NaN;
    lat(lat <= -9999) = NaN;
    lon(lon <= -9999) = NaN;

    % 经度统一到 [-180, 180)
    validLon = isfinite(lon);
    lon(validLon) = mod(lon(validLon) + 180, 360) - 180;

    % 只保留已知天线编码
    ant(~ismember(ant, [2, 3])) = NaN;

    % RFI 标记若不是 0/1，可按需进一步清理
    rfi(~ismember(rfi, [0, 1])) = NaN;

    % Kurtosis 合法性
    kurt(~isfinite(kurt)) = NaN;

    for i = 1:numel(ROIs)
        inROI = lat >= ROIs(i).latMin & lat < ROIs(i).latMax & ...
                lon >= ROIs(i).lonMin & lon < ROIs(i).lonMax;

        % 以 noise + lat + lon + ant 有效作为基础有效样本
        baseValid = inROI & isfinite(noise) & isfinite(lat) & isfinite(lon) & isfinite(ant);

        if ~any(baseValid)
            continue;
        end

        % 逐指标分别统计，避免某一个变量缺测把其他变量也拖掉
        noiseValid = baseValid & isfinite(noise);
        rfiValid   = baseValid & isfinite(rfi);
        kurtValid  = baseValid & isfinite(kurt);

        portMask = baseValid & (ant == 3);   % 3 = Port
        starMask = baseValid & (ant == 2);   % 2 = Starboard

        recIdx = recIdx + 1;
        rec(recIdx).roi   = string(ROIs(i).name);
        rec(recIdx).month = monthDT;

        rec(recIdx).noise_sum = sum(noise(noiseValid), 'omitnan');
        rec(recIdx).noise_n   = sum(noiseValid);

        rec(recIdx).rfi_sum = sum(rfi(rfiValid) == 1);
        rec(recIdx).rfi_n   = sum(rfiValid);

        rec(recIdx).kurt_sum = sum(abs(kurt(kurtValid) - 3), 'omitnan');
        rec(recIdx).kurt_n   = sum(kurtValid);

        rec(recIdx).total_n = sum(baseValid);

        rec(recIdx).port_noise_sum = sum(noise(portMask), 'omitnan');
        rec(recIdx).port_n         = sum(portMask);

        rec(recIdx).star_noise_sum = sum(noise(starMask), 'omitnan');
        rec(recIdx).star_n         = sum(starMask);
    end

    if mod(k, 100) == 0 || k == numel(fileList)
        fprintf('已处理 %d / %d 个文件\n', k, numel(fileList));
    end
end

if isempty(rec)
    error('没有得到任何统计结果。请检查 ROI、变量名、文件日期格式。');
end

%% ===================== 5) 月度汇总 =====================
T = struct2table(rec);

[G, roiG, monthG] = findgroups(T.roi, T.month);

S = table;
S.roi   = splitapply(@(x) x(1), T.roi, G);
S.month = splitapply(@(x) x(1), T.month, G);

S.noise_sum = splitapply(@sum, T.noise_sum, G);
S.noise_n   = splitapply(@sum, T.noise_n, G);

S.rfi_sum = splitapply(@sum, T.rfi_sum, G);
S.rfi_n   = splitapply(@sum, T.rfi_n, G);

S.kurt_sum = splitapply(@sum, T.kurt_sum, G);
S.kurt_n   = splitapply(@sum, T.kurt_n, G);

S.total_n = splitapply(@sum, T.total_n, G);

S.port_noise_sum = splitapply(@sum, T.port_noise_sum, G);
S.port_n         = splitapply(@sum, T.port_n, G);

S.star_noise_sum = splitapply(@sum, T.star_noise_sum, G);
S.star_n         = splitapply(@sum, T.star_n, G);

%% ===================== 6) 计算最终指标 =====================
S.mean_noise_floor = S.noise_sum ./ S.noise_n;
S.rfi_rate         = S.rfi_sum   ./ S.rfi_n;
S.mean_abs_kurt3   = S.kurt_sum  ./ S.kurt_n;

S.port_mean = S.port_noise_sum ./ S.port_n;
S.star_mean = S.star_noise_sum ./ S.star_n;
S.port_minus_star = S.port_mean - S.star_mean;

% 样本太少时，不给 Port-Starboard 结论
badAnt = (S.port_n < minAntCount) | (S.star_n < minAntCount);
S.port_minus_star(badAnt) = NaN;

% Port / Star 样本平衡比，排错很重要
S.port_star_ratio = S.port_n ./ S.star_n;

%% ===================== 7) 排序与导出 =====================
S = sortrows(S, {'roi','month'});
outFile = fullfile(dataDir, 'CYGNSS_ROI_monthly_stats.csv');
writetable(S, outFile);

fprintf('\n月统计表已保存到:\n%s\n', outFile);

disp(S(:, {'roi','month','mean_noise_floor','rfi_rate','mean_abs_kurt3', ...
           'port_minus_star','total_n','port_n','star_n','port_star_ratio'}));

%% ===================== 8) 画图 =====================
roiNames = unique(S.roi);

for i = 1:numel(roiNames)
    idx = S.roi == roiNames(i);
    A = S(idx,:);
    A = sortrows(A, 'month');

    figure('Name', char(roiNames(i)), 'Color', 'w');
    tiledlayout(5,1, 'Padding','compact', 'TileSpacing','compact');

    nexttile;
    plot(A.month, A.mean_noise_floor, '-o', 'LineWidth', 1.2);
    ylabel('mean noise');
    title(char(roiNames(i)));
    grid on;

    nexttile;
    plot(A.month, A.rfi_rate, '-o', 'LineWidth', 1.2);
    ylabel('rfi rate');
    grid on;

    nexttile;
    plot(A.month, A.mean_abs_kurt3, '-o', 'LineWidth', 1.2);
    ylabel('|kurt-3|');
    grid on;

    nexttile;
    plot(A.month, A.port_minus_star, '-o', 'LineWidth', 1.2);
    ylabel('P - S');
    grid on;

    nexttile;
    yyaxis left;
    plot(A.month, A.total_n, '-o', 'LineWidth', 1.2);
    ylabel('total n');

    yyaxis right;
    plot(A.month, A.port_star_ratio, '-s', 'LineWidth', 1.2);
    ylabel('P/S count');
    yline(1, '--');
    grid on;
    xlabel('month');
end



%% ===================== 本地函数 =====================
function dt = parseDateFromFilename(fname)
    % 尝试从文件名中提取 8 位日期，如 20210415
    tk = regexp(fname, '(20\d{6})', 'tokens', 'once');
    if isempty(tk)
        dt = NaT;
        return;
    end
    try
        dt = datetime(tk{1}, 'InputFormat', 'yyyyMMdd');
    catch
        dt = NaT;
    end
ends);

%% ===================== 本地函数 =====================
function dt = parseDateFromFilename(fname)
    % 尝试从文件名中提取 8 位日期，如 20210415
    tk = regexp(fname, '(20\d{6})', 'tokens', 'once');
    if isempty(tk)
        dt = NaT;
        return;
    end
    try
        dt = datetime(tk{1}, 'InputFormat', 'yyyyMMdd');
    catch
        dt = NaT;
    end
endd