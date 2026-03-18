%% ================================================================
%  CYGNSS Port / Starboard 后处理分析脚本
%  用途：
%    1) 计算按纬度的 Port / Starboard 平均曲线
%    2) 计算多个热点 ROI 的统计表
%    3) 输出图和 CSV 文件
%
%  使用方法：
%    先运行 Gleason_fig2_fixed.m，确保以下变量已在工作区中：
%      lat_c, lon_c,
%      nf_mean_sb, nf_mean_pt,
%      cnt_sb, cnt_pt
%    然后再运行本脚本。
%% ================================================================

clc;

%% ---------------- 前置检查 ----------------
required_vars = {'lat_c','lon_c','nf_mean_sb','nf_mean_pt','cnt_sb','cnt_pt'};
for k = 1:numel(required_vars)
    if ~evalin('base', sprintf("exist('%s','var')", required_vars{k}))
        error(['缺少变量: ', required_vars{k}, '。请先运行 Gleason_fig2_fixed.m，再运行本脚本。']);
    end
end

% 从 base workspace 取变量
lat_c      = evalin('base','lat_c');
lon_c      = evalin('base','lon_c');
nf_mean_sb = evalin('base','nf_mean_sb');
nf_mean_pt = evalin('base','nf_mean_pt');
cnt_sb     = evalin('base','cnt_sb');
cnt_pt     = evalin('base','cnt_pt');

%% ---------------- 参数 ----------------
MIN_SAMP_GRID = 30;   % 差值图/ROI分析时，每根天线每格至少样本数
MIN_SAMP_LAT  = 200;  % 按纬度统计时，该纬带每根天线至少总样本数
SAVE_PREFIX   = 'CYGNSS_Port_Starboard';

% 可按你的研究需要修改 ROI 范围
roi_names = {
    'MiddleEast_NorthAfrica'
    'EastAsia_NWPacific'
    'SouthIndianOcean'
    'SouthPacific'
    };

roi_lat_min = [ 20,  20, -35, -35];
roi_lat_max = [ 40,  40, -10, -10];
roi_lon_min = [ 10, 115,  70,-130];
roi_lon_max = [ 60, 150, 110, -95];

%% ---------------- 低样本格点掩膜 ----------------
valid_grid = (cnt_sb >= MIN_SAMP_GRID) & (cnt_pt >= MIN_SAMP_GRID) & ...
             isfinite(nf_mean_sb) & isfinite(nf_mean_pt);

nf_mean_sb_masked = nf_mean_sb;
nf_mean_pt_masked = nf_mean_pt;
diff_map_masked   = nf_mean_pt - nf_mean_sb;

nf_mean_sb_masked(~valid_grid) = NaN;
nf_mean_pt_masked(~valid_grid) = NaN;
diff_map_masked(~valid_grid)   = NaN;

%% ================================================================
%  Part A：按纬度的加权平均曲线
%% ================================================================

n_lat = numel(lat_c);

zonal_sb      = nan(n_lat,1);
zonal_pt      = nan(n_lat,1);
zonal_diff    = nan(n_lat,1);
zonal_cnt_sb  = zeros(n_lat,1);
zonal_cnt_pt  = zeros(n_lat,1);
zonal_ngrid   = zeros(n_lat,1);

for i = 1:n_lat
    sb_row   = nf_mean_sb_masked(i,:);
    pt_row   = nf_mean_pt_masked(i,:);
    sb_cnt   = cnt_sb(i,:);
    pt_cnt   = cnt_pt(i,:);
    row_mask = valid_grid(i,:);

    zonal_ngrid(i) = sum(row_mask);

    if any(row_mask)
        sb_w = sb_cnt(row_mask);
        pt_w = pt_cnt(row_mask);
        sb_v = sb_row(row_mask);
        pt_v = pt_row(row_mask);

        zonal_cnt_sb(i) = sum(sb_w);
        zonal_cnt_pt(i) = sum(pt_w);

        if zonal_cnt_sb(i) >= MIN_SAMP_LAT
            zonal_sb(i) = sum(sb_v .* sb_w) / sum(sb_w);
        end
        if zonal_cnt_pt(i) >= MIN_SAMP_LAT
            zonal_pt(i) = sum(pt_v .* pt_w) / sum(pt_w);
        end
        if isfinite(zonal_sb(i)) && isfinite(zonal_pt(i))
            zonal_diff(i) = zonal_pt(i) - zonal_sb(i);
        end
    end
end

% 保存纬向统计表
T_lat = table(lat_c(:), zonal_sb, zonal_pt, zonal_diff, zonal_cnt_sb, zonal_cnt_pt, zonal_ngrid, ...
    'VariableNames', {'lat_deg','mean_sb_dBW','mean_pt_dBW','pt_minus_sb_dBW', ...
                      'sample_sb','sample_pt','num_valid_grids'});
writetable(T_lat, [SAVE_PREFIX '_latitudinal_profile.csv']);
fprintf('已保存：%s\n', [SAVE_PREFIX '_latitudinal_profile.csv']);

% 图1：Port 和 Starboard 纬向均值
figure('Position',[100,100,860,520]);
plot(lat_c, zonal_sb, 'LineWidth', 2); hold on;
plot(lat_c, zonal_pt, 'LineWidth', 2);
grid on; box on;
xlabel('Latitude (°)');
ylabel('Weighted Mean Noise Floor (dBW)');
legend('Starboard','Port','Location','best');
title('Latitudinal Weighted Mean Noise Floor by Antenna');
saveas(gcf, [SAVE_PREFIX '_latitudinal_mean_sb_pt.png']);
fprintf('已保存：%s\n', [SAVE_PREFIX '_latitudinal_mean_sb_pt.png']);

% 图2：Port - Starboard 纬向差值
figure('Position',[100,100,860,520]);
plot(lat_c, zonal_diff, 'LineWidth', 2); hold on;
yline(0,'--k','LineWidth',1);
grid on; box on;
xlabel('Latitude (°)');
ylabel('Port - Starboard (dBW)');
title('Latitudinal Weighted Mean of Port - Starboard');
saveas(gcf, [SAVE_PREFIX '_latitudinal_diff.png']);
fprintf('已保存：%s\n', [SAVE_PREFIX '_latitudinal_diff.png']);

%% ================================================================
%  Part B：热点 ROI 统计表
%% ================================================================

[LonGrid, LatGrid] = meshgrid(lon_c, lat_c);

roi_mean_sb      = nan(numel(roi_names),1);
roi_mean_pt      = nan(numel(roi_names),1);
roi_mean_diff    = nan(numel(roi_names),1);
roi_num_grids    = zeros(numel(roi_names),1);
roi_sample_sb    = zeros(numel(roi_names),1);
roi_sample_pt    = zeros(numel(roi_names),1);
roi_lat_range    = strings(numel(roi_names),1);
roi_lon_range    = strings(numel(roi_names),1);

for r = 1:numel(roi_names)
    roi_mask = valid_grid & ...
               LatGrid >= roi_lat_min(r) & LatGrid < roi_lat_max(r) & ...
               LonGrid >= roi_lon_min(r) & LonGrid < roi_lon_max(r);

    roi_num_grids(r) = sum(roi_mask(:));
    roi_lat_range(r) = sprintf('[%g, %g)', roi_lat_min(r), roi_lat_max(r));
    roi_lon_range(r) = sprintf('[%g, %g)', roi_lon_min(r), roi_lon_max(r));

    if roi_num_grids(r) == 0
        continue;
    end

    sb_vals = nf_mean_sb_masked(roi_mask);
    pt_vals = nf_mean_pt_masked(roi_mask);
    sb_w    = cnt_sb(roi_mask);
    pt_w    = cnt_pt(roi_mask);

    roi_sample_sb(r) = sum(sb_w);
    roi_sample_pt(r) = sum(pt_w);

    if roi_sample_sb(r) > 0
        roi_mean_sb(r) = sum(sb_vals .* sb_w) / sum(sb_w);
    end
    if roi_sample_pt(r) > 0
        roi_mean_pt(r) = sum(pt_vals .* pt_w) / sum(pt_w);
    end
    if isfinite(roi_mean_sb(r)) && isfinite(roi_mean_pt(r))
        roi_mean_diff(r) = roi_mean_pt(r) - roi_mean_sb(r);
    end
end

T_roi = table(roi_names, roi_lat_range, roi_lon_range, roi_num_grids, ...
              roi_sample_sb, roi_sample_pt, ...
              roi_mean_sb, roi_mean_pt, roi_mean_diff, ...
    'VariableNames', {'ROI','lat_range_deg','lon_range_deg','num_valid_grids', ...
                      'sample_sb','sample_pt', ...
                      'mean_sb_dBW','mean_pt_dBW','pt_minus_sb_dBW'});

writetable(T_roi, [SAVE_PREFIX '_ROI_stats.csv']);
disp(T_roi);
fprintf('已保存：%s\n', [SAVE_PREFIX '_ROI_stats.csv']);

%% ================================================================
%  Part C：可视化 ROI 位置
%% ================================================================
load coastlines;
figure('Position',[80,80,1200,450]);
imagesc(lon_c, lat_c, diff_map_masked);
set(gca,'YDir','normal');
hold on;
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
colormap(redblue_colormap());
cb = colorbar;
cb.Label.String = 'Port - Starboard (dBW)';
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('ROI Boxes on Port - Starboard Difference Map');
grid on; box on;

for r = 1:numel(roi_names)
    rectangle('Position', [roi_lon_min(r), roi_lat_min(r), ...
                           roi_lon_max(r)-roi_lon_min(r), ...
                           roi_lat_max(r)-roi_lat_min(r)], ...
              'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');
    text(roi_lon_min(r)+1, roi_lat_max(r)-1, roi_names{r}, ...
         'Color','k','FontSize',9,'FontWeight','bold', ...
         'BackgroundColor','w','Margin',1);
end
saveas(gcf, [SAVE_PREFIX '_ROI_boxes_on_diffmap.png']);
fprintf('已保存：%s\n', [SAVE_PREFIX '_ROI_boxes_on_diffmap.png']);

fprintf('\n后处理完成。\n');

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