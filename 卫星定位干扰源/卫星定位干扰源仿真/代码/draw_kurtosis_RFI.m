%% ================================================
%  Step 1: RFI自动识别与特征关联
%  基于 ddm_kurtosis 的 RFI 检测，对标 Al-Khaldi 2025
%  在原有 Figure 1 框架上扩展
%% ================================================
clear; clc; close all;

%% ---- 参数（与原代码保持一致）----
base_dir  = 'D:\CYGNSS_L1_v32\2021\2021.01\';
sat_list  = {'cyg01','cyg02','cyg03','cyg05','cyg06','cyg07','cyg08'};
RFI_MASK  = uint32(131072);   % quality_flags 第18位 = rfi_detected
KAPPA_THR = 1.0;              % 官方L1产品阈值：Δκ > 1.0 判为RFI（与quality_flag一致）
                              % 注：Al-Khaldi 2025的0.1阈值针对原始I/F数据，非此处L1产品
lat_edges = -40:1:40;
lon_edges = -180:1:180;
n_lat     = length(lat_edges) - 1;
n_lon     = length(lon_edges) - 1;
N         = n_lat * n_lon;
MIN_SAMP  = 10;

%% ---- 初始化累加器 ----
nf_sum       = zeros(n_lat, n_lon);   % noise floor 累加（dB）
dk_sum       = zeros(n_lat, n_lon);   % Δκ 累加
rfi_flag_sum = zeros(n_lat, n_lon);   % quality_flag RFI 计数
rfi_kur_sum  = zeros(n_lat, n_lon);   % kurtosis RFI 计数
cnt_total    = zeros(n_lat, n_lon);   % 总采样计数

% 用于散点图的随机采样池（避免内存爆炸）
MAX_SCATTER  = 200000;
scatter_nf   = nan(MAX_SCATTER, 1);
scatter_dk   = nan(MAX_SCATTER, 1);
scatter_flag = false(MAX_SCATTER, 1);
scatter_ptr  = 0;

%% ---- 遍历所有日期文件夹 ----
day_folders = dir(fullfile(base_dir, '2021.01.*'));
fprintf('找到 %d 个日期文件夹\n', length(day_folders));

for d = 1:length(day_folders)
    day_dir = fullfile(base_dir, day_folders(d).name);
    fprintf('\n[%d/%d] %s\n', d, length(day_folders), day_folders(d).name);

    for s = 1:length(sat_list)
        files = dir(fullfile(day_dir, [sat_list{s}, '*.nc']));
        if isempty(files), continue; end
        filename = fullfile(day_dir, files(1).name);

        % ---- 读取变量（新增 ddm_kurtosis）----
        try
            sp_lat      = ncread(filename, 'sp_lat');
            sp_lon      = ncread(filename, 'sp_lon');
            noise_floor = ncread(filename, 'ddm_noise_floor');
            qflags      = ncread(filename, 'quality_flags');
            kurtosis_raw = ncread(filename, 'ddm_kurtosis');  % ★ 新增
        catch ME
            fprintf('  读取失败：%s\n  原因：%s\n', files(1).name, ME.message);
            continue
        end

        % ---- 清洗 ----
        noise_floor(noise_floor <= 0)    = NaN;
        kurtosis_raw(kurtosis_raw <= -9) = NaN;  % 过滤填充值-9999及其他异常值
        sp_lat(sp_lat  <= -9999)        = NaN;
        sp_lon(sp_lon  <= -9999)        = NaN;
        sp_lon(sp_lon > 180)            = sp_lon(sp_lon > 180) - 360;

        % ---- 计算 Δκ = |κ - 3| ----
        % 理论高斯白噪声峰度 = 3，偏离越大越可能是RFI
        delta_kappa = abs(kurtosis_raw - 3);

        % ---- 提取标志 ----
        rfi_flag    = bitand(uint32(qflags), RFI_MASK) > 0;   % 官方标记
        rfi_kur_det = delta_kappa > KAPPA_THR;                 % 峰度检测

        % ---- 展平为向量 ----
        lat_vec  = sp_lat(:);
        lon_vec  = sp_lon(:);
        nf_vec   = 10 * log10(noise_floor(:));
        dk_vec   = delta_kappa(:);
        flag_vec = double(rfi_flag(:));
        kur_vec  = double(rfi_kur_det(:));

        valid = ~isnan(lat_vec) & ~isnan(lon_vec) & ...
                ~isnan(nf_vec)  & ~isnan(dk_vec);

        lat_v  = lat_vec(valid);
        lon_v  = lon_vec(valid);
        nf_v   = nf_vec(valid);
        dk_v   = dk_vec(valid);
        flag_v = flag_vec(valid);
        kur_v  = kur_vec(valid);

        % ---- 散点图随机采样 ----
        n_new = sum(valid);
        if n_new > 0 && scatter_ptr < MAX_SCATTER
            idx_sample = randperm(n_new, min(n_new, MAX_SCATTER - scatter_ptr));
            n_add = length(idx_sample);
            scatter_nf  (scatter_ptr+1 : scatter_ptr+n_add) = nf_v(idx_sample);
            scatter_dk  (scatter_ptr+1 : scatter_ptr+n_add) = dk_v(idx_sample);
            scatter_flag(scatter_ptr+1 : scatter_ptr+n_add) = logical(flag_v(idx_sample));
            scatter_ptr = scatter_ptr + n_add;
        end

        % ---- 空间格网累加 ----
        [~, lat_idx] = histc(lat_v, lat_edges);
        [~, lon_idx] = histc(lon_v, lon_edges);
        ok = lat_idx>0 & lat_idx<=n_lat & lon_idx>0 & lon_idx<=n_lon;
        if sum(ok) == 0, continue; end

        idx = sub2ind([n_lat, n_lon], lat_idx(ok), lon_idx(ok));

        nf_sum       = nf_sum       + reshape(accumarray(idx, nf_v(ok),   [N,1],@sum,0), [n_lat,n_lon]);
        dk_sum       = dk_sum       + reshape(accumarray(idx, dk_v(ok),   [N,1],@sum,0), [n_lat,n_lon]);
        rfi_flag_sum = rfi_flag_sum + reshape(accumarray(idx, flag_v(ok), [N,1],@sum,0), [n_lat,n_lon]);
        rfi_kur_sum  = rfi_kur_sum  + reshape(accumarray(idx, kur_v(ok),  [N,1],@sum,0), [n_lat,n_lon]);
        cnt_total    = cnt_total    + reshape(accumarray(idx, ones(sum(ok),1), [N,1],@sum,0), [n_lat,n_lon]);
    end

    fprintf('  累计有效采样：%.2fM\n', sum(cnt_total(:))/1e6);
end

%% ---- 计算均值 ----
nf_mean      = nf_sum       ./ max(cnt_total, 1);
dk_mean      = dk_sum       ./ max(cnt_total, 1);
rfi_flag_rate = rfi_flag_sum ./ max(cnt_total, 1) * 100;
rfi_kur_rate  = rfi_kur_sum  ./ max(cnt_total, 1) * 100;

% 采样不足的格子置NaN
mask = cnt_total < MIN_SAMP;
nf_mean(mask)       = NaN;
dk_mean(mask)       = NaN;
rfi_flag_rate(mask) = NaN;
rfi_kur_rate(mask)  = NaN;

lat_centers = lat_edges(1:end-1) + 0.5;
lon_centers = lon_edges(1:end-1) + 0.5;

%% ---- 统计输出 ----
fprintf('\n==============================\n');
fprintf('官方RFI标记率：%.2f%%\n',  100*sum(rfi_flag_sum(:))/sum(cnt_total(:)));
fprintf('峰度RFI检测率：%.2f%%\n',  100*sum(rfi_kur_sum(:)) /sum(cnt_total(:)));
fprintf('总采样点数：%.2fM\n', sum(cnt_total(:))/1e6);

%% ---- 两种检测方法混淆矩阵 ----
% 利用散点数据（逐点对比）
valid_sc = ~isnan(scatter_nf(1:scatter_ptr)) & ~isnan(scatter_dk(1:scatter_ptr));
sc_flag  = scatter_flag(1:scatter_ptr);
sc_kur   = scatter_dk(1:scatter_ptr) > KAPPA_THR;
sc_valid_flag = sc_flag(valid_sc);
sc_valid_kur  = sc_kur(valid_sc);

TP = sum( sc_valid_flag &  sc_valid_kur);
FP = sum(~sc_valid_flag &  sc_valid_kur);
FN = sum( sc_valid_flag & ~sc_valid_kur);
TN = sum(~sc_valid_flag & ~sc_valid_kur);
fprintf('\n--- 峰度检测 vs 官方标记 混淆矩阵（采样 %.0fK 点）---\n', scatter_ptr/1000);
fprintf('              官方=RFI   官方=正常\n');
fprintf('峰度检测=RFI  %8d   %8d\n', TP, FP);
fprintf('峰度检测=正常 %8d   %8d\n', FN, TN);
prec = TP/(TP+FP+eps); rec = TP/(TP+FN+eps);
fprintf('精确率 Precision: %.1f%%   召回率 Recall: %.1f%%\n', prec*100, rec*100);

%% ================================================================
%  图1：全局对比图（4子图）
%% ================================================================
load coastlines;
figure('Position',[30,30,1500,900], 'Name','全球RFI分析');

% -- 子图1：噪声底均值 --
subplot(2,2,1);
imagesc(lon_centers, lat_centers, nf_mean);
set(gca,'YDir','normal'); colorbar; colormap(gca, jet(256));
nf_v2 = nf_mean(~isnan(nf_mean));
clim([prctile(nf_v2,2), prctile(nf_v2,98)]);
hold on; plot(coastlon,coastlat,'k-','LineWidth',0.5); hold off;
title('(a) 平均噪声底 (dB)');
xlabel('Longitude (°)'); ylabel('Latitude (°)');

% -- 子图2：官方RFI标记率 --
subplot(2,2,2);
imagesc(lon_centers, lat_centers, rfi_flag_rate);
set(gca,'YDir','normal'); colorbar; colormap(gca, hot(256));
clim([0, min(100, prctile(rfi_flag_rate(~isnan(rfi_flag_rate)),99))]);
hold on; plot(coastlon,coastlat,'w-','LineWidth',0.5);
plot(35.95,35.4,'c^','MarkerSize',9,'MarkerFaceColor','c');
plot(114,32,'g^','MarkerSize',9,'MarkerFaceColor','g');
text(37,34,'Syria','Color','cyan','FontSize',8,'FontWeight','bold');
text(116,31,'China','Color','green','FontSize',8,'FontWeight','bold');
hold off;
title('(b) 官方 quality\_flag RFI率 (%)');
xlabel('Longitude (°)'); ylabel('Latitude (°)');

% -- 子图3：峰度偏离均值图 --
subplot(2,2,3);
imagesc(lon_centers, lat_centers, dk_mean);
set(gca,'YDir','normal'); colorbar; colormap(gca, jet(256));
dk_v2 = dk_mean(~isnan(dk_mean));
clim([0, prctile(dk_v2,98)]);
hold on; plot(coastlon,coastlat,'k-','LineWidth',0.5); hold off;
title(sprintf('(c) 平均 Δκ = |κ−3| (阈值=%.1f)', KAPPA_THR));
xlabel('Longitude (°)'); ylabel('Latitude (°)');

% -- 子图4：峰度检测RFI率 --
subplot(2,2,4);
imagesc(lon_centers, lat_centers, rfi_kur_rate);
set(gca,'YDir','normal'); colorbar; colormap(gca, hot(256));
clim([0, min(100, prctile(rfi_kur_rate(~isnan(rfi_kur_rate)),99))]);
hold on; plot(coastlon,coastlat,'w-','LineWidth',0.5);
plot(35.95,35.4,'c^','MarkerSize',9,'MarkerFaceColor','c');
plot(114,32,'g^','MarkerSize',9,'MarkerFaceColor','g');
text(37,34,'Syria','Color','cyan','FontSize',8,'FontWeight','bold');
text(116,31,'China','Color','green','FontSize',8,'FontWeight','bold');
hold off;
title(sprintf('(d) 峰度检测 RFI率 (%%) [Δκ>%.1f]', KAPPA_THR));
xlabel('Longitude (°)'); ylabel('Latitude (°)');

sgtitle('CYGNSS L1 v3.2 — 2021年1月 RFI检测对比（官方标记 vs 峰度检测）', ...
        'FontSize',13,'FontWeight','bold');
saveas(gcf, 'CYGNSS_kurtosis_vs_flag_global.png');

%% ================================================================
%  图2：散点图 — Δκ vs Noise Floor，按官方RFI标记着色
%% ================================================================
valid_sc = ~isnan(scatter_nf(1:scatter_ptr)) & ~isnan(scatter_dk(1:scatter_ptr));
sc_nf    = scatter_nf(1:scatter_ptr);
sc_dk    = scatter_dk(1:scatter_ptr);
sc_f     = scatter_flag(1:scatter_ptr);

% 限制散点数量，提高绘图速度
MAX_PLOT  = 50000;
idx_plot  = find(valid_sc);
if length(idx_plot) > MAX_PLOT
    idx_plot = idx_plot(randperm(length(idx_plot), MAX_PLOT));
end

figure('Position',[100,100,900,600], 'Name','散点图：Δκ vs 噪声底');
hold on;
% 官方标记为正常的点（蓝色，小）
idx_norm = idx_plot(~sc_f(idx_plot));
idx_rfi  = idx_plot( sc_f(idx_plot));
scatter(sc_nf(idx_norm), sc_dk(idx_norm), 4, [0.4 0.6 1.0], 'filled', 'MarkerFaceAlpha', 0.3);
scatter(sc_nf(idx_rfi),  sc_dk(idx_rfi),  8, [1.0 0.3 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
% 峰度阈值线
yline(KAPPA_THR, 'k--', 'LineWidth',1.5, ...
      'Label', sprintf('Δκ = %.1f 阈值', KAPPA_THR), 'LabelHorizontalAlignment','left');
hold off;
xlabel('噪声底 (dB)');
ylabel('Δκ = |κ − 3|');
legend('官方标记：正常','官方标记：RFI','Location','northeast');
title('Δκ vs 噪声底 散点图（2021年1月随机采样）');
set(gca,'YScale','log');
grid on;
saveas(gcf, 'CYGNSS_scatter_dk_vs_nf.png');

%% ================================================================
%  图3：Δκ 累积分布函数（对标 Al-Khaldi 2025 Fig.5）
%% ================================================================
dk_all = sc_dk(valid_sc);
dk_all = dk_all(dk_all >= 0 & dk_all < 10);  % 去除极端异常值

thr_vec   = linspace(0, 2.0, 400);
ccdf_vals = arrayfun(@(t) sum(dk_all > t)/length(dk_all)*100, thr_vec);

figure('Position',[150,150,700,500], 'Name','Δκ 互补累积分布');
plot(thr_vec, ccdf_vals, 'b-', 'LineWidth',2);
hold on;
xline(KAPPA_THR, 'r--', 'LineWidth',1.5, ...
      'Label', sprintf('%.1f (%.1f%%)', KAPPA_THR, ...
      interp1(thr_vec,ccdf_vals,KAPPA_THR,'linear')), ...
      'LabelHorizontalAlignment','right');
hold off;
xlabel('Δκ 阈值');
ylabel('检测比例 (%)');
title('互补累积分布函数 1-CDF（对标 Al-Khaldi 2025 Fig.5）');
grid on; ylim([0, 100]); xlim([0, 2.0]);
saveas(gcf, 'CYGNSS_kurtosis_CCDF.png');

fprintf('\n所有图像已保存。\n');
fprintf('建议对比：\n');
fprintf('  - 图(b) vs 图(d)：官方标记 与 峰度检测 的空间一致性\n');
fprintf('  - 散点图中红点（官方RFI）是否集中在高Δκ区域\n');
fprintf('  - CCDF图与 Al-Khaldi 2025 Fig.5 对比阈值选取\n');
