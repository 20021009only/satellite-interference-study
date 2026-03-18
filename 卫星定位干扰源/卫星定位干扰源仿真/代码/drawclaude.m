%% ================================================
%  多天多星合并：月平均全球RFI图
%  对标论文 Fig.1
%% ================================================
clear; clc; close all;

%% ---- 参数 ----
base_dir  = 'D:\CYGNSS_L1_v32\2021\2021.01\';
sat_list  = {'cyg01','cyg02','cyg03','cyg05','cyg06','cyg07','cyg08'};
RFI_MASK  = uint32(131072);
lat_edges = -40:1:40;
lon_edges = -180:1:180;
n_lat     = length(lat_edges)-1;
n_lon     = length(lon_edges)-1;
N         = n_lat * n_lon;
MIN_SAMP  = 10;

%% ---- 初始化 ----
nf_sum    = zeros(n_lat, n_lon);
rfi_sum   = zeros(n_lat, n_lon);
cnt_total = zeros(n_lat, n_lon);

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

        try
            sp_lat      = ncread(filename, 'sp_lat');
            sp_lon      = ncread(filename, 'sp_lon');
            noise_floor = ncread(filename, 'ddm_noise_floor');
            qflags      = ncread(filename, 'quality_flags');
        catch
            fprintf('  读取失败：%s\n', files(1).name);
            continue
        end

        noise_floor(noise_floor <= 0) = NaN;
        sp_lat(sp_lat   <= -9999)     = NaN;
        sp_lon(sp_lon   <= -9999)     = NaN;
        sp_lon(sp_lon > 180)          = sp_lon(sp_lon > 180) - 360;

        rfi_flag = bitand(uint32(qflags), RFI_MASK) > 0;

        lat_vec = sp_lat(:);
        lon_vec = sp_lon(:);
        nf_vec  = noise_floor(:);
        rfi_vec = double(rfi_flag(:));

        valid = ~isnan(lat_vec) & ~isnan(lon_vec) & ~isnan(nf_vec);
        lat_vec = lat_vec(valid);
        lon_vec = lon_vec(valid);
        nf_vec  = 10*log10(nf_vec(valid));
        rfi_vec = rfi_vec(valid);

        [~, lat_idx] = histc(lat_vec, lat_edges);
        [~, lon_idx] = histc(lon_vec, lon_edges);
        ok = lat_idx>0 & lat_idx<=n_lat & lon_idx>0 & lon_idx<=n_lon;
        if sum(ok) == 0, continue; end

        idx = sub2ind([n_lat,n_lon], lat_idx(ok), lon_idx(ok));

        nf_sum    = nf_sum    + reshape(accumarray(idx,nf_vec(ok),       [N,1],@sum,0),[n_lat,n_lon]);
        rfi_sum   = rfi_sum   + reshape(accumarray(idx,rfi_vec(ok),      [N,1],@sum,0),[n_lat,n_lon]);
        cnt_total = cnt_total + reshape(accumarray(idx,ones(sum(ok),1),  [N,1],@sum,0),[n_lat,n_lon]);
    end

    % 进度显示
    fprintf('  累计有效采样：%.2fM\n', sum(cnt_total(:))/1e6);
end

%% ---- 计算均值并过滤 ----
nf_mean  = nf_sum  ./ max(cnt_total, 1);
rfi_rate = rfi_sum ./ max(cnt_total, 1) * 100;
nf_mean(cnt_total  < MIN_SAMP) = NaN;
rfi_rate(cnt_total < MIN_SAMP) = NaN;

lat_centers = lat_edges(1:end-1) + 0.5;
lon_centers = lon_edges(1:end-1) + 0.5;

fprintf('\n==============================\n');
fprintf('全局平均RFI率：%.2f%%\n', 100*sum(rfi_sum(:))/sum(cnt_total(:)));
fprintf('总采样点数：%.2fM\n', sum(cnt_total(:))/1e6);

%% ---- 画图 ----
figure('Position',[50,50,1400,800]);

subplot(2,1,1);
imagesc(lon_centers, lat_centers, nf_mean);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'Mean Noise Floor (dB)';
colormap(gca, jet(256));
xlabel('经度 (°E)'); ylabel('纬度 (°N)');
title('CYGNSS L1 v3.2 全球噪声底本均值 — 2021年1月（7星×31天）');
nf_vals = nf_mean(~isnan(nf_mean));
clim([prctile(nf_vals,2), prctile(nf_vals,98)]);
hold on; load coastlines;
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.5); hold off;

subplot(2,1,2);
imagesc(lon_centers, lat_centers, rfi_rate);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'RFI Detection Rate (%)';
colormap(gca, hot(256));
xlabel('经度 (°E)'); ylabel('纬度 (°N)');
title('CYGNSS L1 v3.2 RFI发生率 — 2021年1月（7星×31天）');
clim([0, 100]);
hold on; load coastlines;
plot(coastlon, coastlat, 'w-', 'LineWidth', 0.5);
% 标注已知干扰源
plot(35.95, 35.4,  'c^', 'MarkerSize',10, 'MarkerFaceColor','c');
plot(114,   32,    'g^', 'MarkerSize',10, 'MarkerFaceColor','g');
text(37, 34,   '叙利亚', 'Color','cyan',  'FontSize',9, 'FontWeight','bold');
text(116, 31,  '中国',   'Color','green', 'FontSize',9, 'FontWeight','bold');
hold off;

saveas(gcf, 'CYGNSS_global_RFI_2021Jan.png');
fprintf('图像已保存\n');