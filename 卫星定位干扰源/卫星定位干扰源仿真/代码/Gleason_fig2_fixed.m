%% ================================================================
%  Gleason 2020 Fig.2 复现（修正版）
%  修复内容：
%  1) 修正填充值过滤 bug：
%     - ddm_noise_floor 有效值应 > 0
%     - sp_lat / sp_lon 只应过滤 _FillValue = -9999
%  2) 保留 ddm_ant 正确编码：2 = Starboard, 3 = Port
%  3) 保留 NaN 灰底显示与样本密度诊断图
%% ================================================================
clear; clc; close all;

%% ---------------- 参数 ----------------
base_dir  = 'D:\CYGNSS_L1_v32\2021\2021.01\';
sat_list  = {'cyg01','cyg02','cyg03','cyg05','cyg06','cyg07','cyg08'};

lat_edges = -40:1:40;
lon_edges = -180:1:180;
n_lat     = length(lat_edges) - 1;
n_lon     = length(lon_edges) - 1;
N         = n_lat * n_lon;
MIN_SAMP  = 10;

% CYGNSS L1 实测编码
ANT_STARBOARD = 2;
ANT_PORT      = 3;

GRAY = [0.75 0.75 0.75];

%% ---------------- 初始化 ----------------
nf_sum_sb = zeros(n_lat, n_lon);
cnt_sb    = zeros(n_lat, n_lon);
nf_sum_pt = zeros(n_lat, n_lon);
cnt_pt    = zeros(n_lat, n_lon);

%% ---------------- 遍历日期文件夹 ----------------
day_folders = dir(fullfile(base_dir, '2021.01.*'));
fprintf('找到 %d 个日期文件夹\n', length(day_folders));

for d = 1:length(day_folders)
    day_dir = fullfile(base_dir, day_folders(d).name);
    fprintf('[%d/%d] %s\n', d, length(day_folders), day_folders(d).name);

    for s = 1:length(sat_list)
        files = dir(fullfile(day_dir, [sat_list{s}, '*.nc']));
        if isempty(files)
            continue;
        end
        filename = fullfile(day_dir, files(1).name);

        try
            sp_lat      = ncread(filename, 'sp_lat');
            sp_lon      = ncread(filename, 'sp_lon');
            noise_floor = ncread(filename, 'ddm_noise_floor');
            ddm_ant     = ncread(filename, 'ddm_ant');
        catch ME
            fprintf('  读取失败：%s (%s)\n', files(1).name, ME.message);
            continue;
        end

        %% ---------------- 关键修复：填充值过滤 ----------------
        % 正确逻辑：
        % 1) ddm_noise_floor 的无效值包含 -9999，且噪声底应为正值，因此 <= 0 都置 NaN
        % 2) sp_lat / sp_lon 的填充值是 -9999，不能写成 <= -9，否则会把所有南纬 9°以南数据误删
        noise_floor(noise_floor <= 0) = NaN;
        sp_lat(sp_lat <= -9999)       = NaN;
        sp_lon(sp_lon <= -9999)       = NaN;

        %% ---------------- 经度转换 ----------------
        % CYGNSS 原始经度为 0~360，需要转为 -180~180
        sp_lon(sp_lon > 180) = sp_lon(sp_lon > 180) - 360;

        %% ---------------- 转 dB ----------------
        % 这里保持你的原统计口径：先转 dB，再做格网平均
        nf_db = 10 .* log10(noise_floor);

        %% ---------------- 展平 ----------------
        lat_v = sp_lat(:);
        lon_v = sp_lon(:);
        nf_v  = nf_db(:);
        ant_v = double(ddm_ant(:));

        %% ---------------- 分天线累计 ----------------
        for ant_code = [ANT_STARBOARD, ANT_PORT]
            mask = (ant_v == ant_code) & ...
                   ~isnan(lat_v) & ~isnan(lon_v) & ...
                   ~isnan(nf_v)  & isfinite(nf_v);

            if ~any(mask)
                continue;
            end

            lv = lat_v(mask);
            lo = lon_v(mask);
            nv = nf_v(mask);

            [~, li] = histc(lv, lat_edges);
            [~, oi] = histc(lo, lon_edges);

            ok = li > 0 & li <= n_lat & oi > 0 & oi <= n_lon;
            if ~any(ok)
                continue;
            end

            idx     = sub2ind([n_lat, n_lon], li(ok), oi(ok));
            nf_acc  = reshape(accumarray(idx, nv(ok),          [N, 1], @sum, 0), [n_lat, n_lon]);
            cnt_acc = reshape(accumarray(idx, ones(sum(ok),1), [N, 1], @sum, 0), [n_lat, n_lon]);

            if ant_code == ANT_STARBOARD
                nf_sum_sb = nf_sum_sb + nf_acc;
                cnt_sb    = cnt_sb    + cnt_acc;
            else
                nf_sum_pt = nf_sum_pt + nf_acc;
                cnt_pt    = cnt_pt    + cnt_acc;
            end
        end
    end
end

%% ---------------- 计算均值 ----------------
nf_mean_sb = nf_sum_sb ./ max(cnt_sb, 1);
nf_mean_pt = nf_sum_pt ./ max(cnt_pt, 1);

nf_mean_sb(cnt_sb < MIN_SAMP) = NaN;
nf_mean_pt(cnt_pt < MIN_SAMP) = NaN;

lat_c = lat_edges(1:end-1) + 0.5;
lon_c = lon_edges(1:end-1) + 0.5;

fprintf('\n==============================\n');
fprintf('Starboard 有效格点：%d，总采样：%.2fM\n', sum(~isnan(nf_mean_sb(:))), sum(cnt_sb(:))/1e6);
fprintf('Port      有效格点：%d，总采样：%.2fM\n', sum(~isnan(nf_mean_pt(:))), sum(cnt_pt(:))/1e6);

%% ---------------- 统一色标 ----------------
all_vals = [nf_mean_sb(~isnan(nf_mean_sb)); nf_mean_pt(~isnan(nf_mean_pt))];
clim_lo  = prctile(all_vals, 2);
clim_hi  = prctile(all_vals, 98);

load coastlines;

%% ================================================================
%  Figure 1：主图（上 Starboard，下 Port）
%% ================================================================
figure('Position', [50, 50, 1200, 720]);

ax1 = subplot(2,1,1);
h1  = imagesc(lon_c, lat_c, nf_mean_sb);
set(h1, 'AlphaData', ~isnan(nf_mean_sb));
set(ax1, 'YDir', 'normal', 'Color', GRAY);
colormap(ax1, jet(256));
clim(ax1, [clim_lo, clim_hi]);
colorbar; xlabel('Longitude (°)'); ylabel('Latitude (°)');
title('(Top) Starboard Mean Noise Floor Power — CYGNSS L1 v3.2, Jan 2021');
hold on; plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

ax2 = subplot(2,1,2);
h2  = imagesc(lon_c, lat_c, nf_mean_pt);
set(h2, 'AlphaData', ~isnan(nf_mean_pt));
set(ax2, 'YDir', 'normal', 'Color', GRAY);
colormap(ax2, jet(256));
clim(ax2, [clim_lo, clim_hi]);
colorbar; xlabel('Longitude (°)'); ylabel('Latitude (°)');
title('(Bottom) Port Mean Noise Floor Power — CYGNSS L1 v3.2, Jan 2021');
hold on; plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

sgtitle('Replication of Gleason 2020 Fig.2: Global DDM Noise Floor by Antenna', ...
        'FontSize', 13, 'FontWeight', 'bold');

saveas(gcf, 'Gleason_Fig2_noise_floor_port_starboard_FIXED.png');
fprintf('主图已保存：Gleason_Fig2_noise_floor_port_starboard_FIXED.png\n');

%% ================================================================
%  Figure 2：差值图（Port - Starboard）
%% ================================================================
diff_map = nf_mean_pt - nf_mean_sb;
diff_vals = diff_map(~isnan(diff_map));
clim_diff = prctile(abs(diff_vals), 98);

figure('Position', [50, 50, 1200, 420]);
ax3 = gca;
h3  = imagesc(lon_c, lat_c, diff_map);
set(h3, 'AlphaData', ~isnan(diff_map));
set(ax3, 'YDir', 'normal', 'Color', GRAY);
colormap(ax3, redblue_colormap());
clim(ax3, [-clim_diff, clim_diff]);
colorbar; xlabel('Longitude (°)'); ylabel('Latitude (°)');
title('Port - Starboard Noise Floor Difference (dBW) — Jan 2021');
hold on; plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

saveas(gcf, 'Gleason_Fig2_diff_port_minus_starboard_FIXED.png');
fprintf('差值图已保存：Gleason_Fig2_diff_port_minus_starboard_FIXED.png\n');

%% ================================================================
%  Figure 3：样本密度图（诊断用）
%% ================================================================
figure('Position', [50, 50, 1200, 720]);

ax4 = subplot(2,1,1);
h4  = imagesc(lon_c, lat_c, log10(cnt_sb + 1));
set(h4, 'AlphaData', cnt_sb > 0);
set(ax4, 'YDir', 'normal', 'Color', GRAY);
colormap(ax4, parula(256));
colorbar; xlabel('Longitude (°)'); ylabel('Latitude (°)');
title('Starboard Sample Density (log scale) — Jan 2021');
hold on; plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

ax5 = subplot(2,1,2);
h5  = imagesc(lon_c, lat_c, log10(cnt_pt + 1));
set(h5, 'AlphaData', cnt_pt > 0);
set(ax5, 'YDir', 'normal', 'Color', GRAY);
colormap(ax5, parula(256));
colorbar; xlabel('Longitude (°)'); ylabel('Latitude (°)');
title('Port Sample Density (log scale) — Jan 2021');
hold on; plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

sgtitle('Sample Density by Antenna — CYGNSS L1 v3.2, Jan 2021', ...
        'FontSize', 13, 'FontWeight', 'bold');

saveas(gcf, 'Gleason_Fig2_sample_density_FIXED.png');
fprintf('样本密度图已保存：Gleason_Fig2_sample_density_FIXED.png\n');

fprintf('\n全部完成。\n');

%% ================================================================
%  辅助函数：红-白-蓝 colormap
%% ================================================================
function cmap = redblue_colormap()
    n = 128;
    r = [linspace(0,1,n), ones(1,n)];
    g = [linspace(0,1,n), linspace(1,0,n)];
    b = [ones(1,n),       linspace(1,0,n)];
    cmap = [r', g', b'];
end
