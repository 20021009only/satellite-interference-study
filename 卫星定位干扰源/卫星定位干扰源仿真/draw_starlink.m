% ========================================================
% 脚本功能：基于 TLE 数据计算星下点轨迹 (Ground Track)
% 作者：Only (Assisted by Gemini)
% ========================================================
clc; clear; close all;

%% 1. 输入 TLE 数据 (就是你刚才读取到的那一行)
% STARLINK-1008 (Line 2)
tle_line2 = '2 44714  53.1593 339.5849 0001299  80.4695 279.6453 15.31021178343579';

% 解析六根数 (单位转换: 角度 -> 弧度)
data = sscanf(tle_line2, '%d %d %f %f %f %f %f %f');
inc = deg2rad(data(3));      % 轨道倾角 (i)
RAAN = deg2rad(data(4));     % 升交点赤经 (Omega)
ecc = data(5) * 1e-7;        % 偏心率 (e) - 注意 TLE 省略了小数点
omega = deg2rad(data(6));    % 近地点幅角 (w)
M0 = deg2rad(data(7));       % 平近点角 (M)
n_rev_per_day = data(8);     % 平均运动 (n) [rev/day]

%% 2. 定义物理常数
mu = 3.986004418e14;         % 地球引力常数 (m^3/s^2)
We = 7.292115e-5;            % 地球自转角速度 (rad/s)
T_day = 86400;               % 一天秒数

% 计算半长轴 a (根据开普勒第三定律: n^2 * a^3 = mu)
n_rad_s = n_rev_per_day * 2 * pi / T_day; % 转换为 rad/s
a = (mu / n_rad_s^2)^(1/3);  % 轨道半长轴 (m)

%% 3. 轨道推演 (Propagate)
% 仿真时间：跑两圈 (约 11200 秒)
t = 0:10:11200; 
lat_store = [];
lon_store = [];

fprintf('正在计算轨道...\n');

for k = 1:length(t)
    dt = t(k);
    
    % A. 更新平近点角 M(t)
    M = M0 + n_rad_s * dt;
    
    % B. 求解开普勒方程 E - e*sin(E) = M 得到偏近点角 E (牛顿迭代法)
    E = M; 
    for iter = 1:10
        E = E - (E - ecc*sin(E) - M) / (1 - ecc*cos(E));
    end
    
    % C. 计算真近点角 (True Anomaly) v
    sin_v = (sqrt(1-ecc^2) * sin(E)) / (1 - ecc*cos(E));
    cos_v = (cos(E) - ecc) / (1 - ecc*cos(E));
    v = atan2(sin_v, cos_v);
    
    % D. 计算轨道平面内的位置 (r) 和 纬度参数 (u)
    r = a * (1 - ecc^2) / (1 + ecc*cos(v));
    u = omega + v; % 纬度幅角
    
    % E. 计算在地心惯性系 (ECI) 中的坐标 (简化版)
    % 这里主要为了求经纬度，不需要完整的 xyz
    % 卫星在轨道面上的投影点纬度 phi
    sin_phi = sin(inc) * sin(u);
    phi = asin(sin_phi); % 纬度
    
    % 卫星在轨道面上的投影点经度 lambda (相对于升交点)
    y_node = cos(inc) * sin(u);
    x_node = cos(u);
    lambda_node = atan2(y_node, x_node);
    
    % F. 转换到地心地固系 (ECEF) -> 也就是经度
    % 核心：必须减去地球自转带来的经度漂移 (We * dt)
    lambda = RAAN + lambda_node - We * dt;
    
    % 存结果 (转为角度)
    lat_store(k) = rad2deg(phi);
    % 处理经度跨越 -180/180 的情况
    lon_deg = mod(rad2deg(lambda), 360);
    if lon_deg > 180, lon_deg = lon_deg - 360; end
    lon_store(k) = lon_deg;
end

%% 4. 绘图：星下点轨迹
figure('Color','w');
plot(lon_store, lat_store, 'b.', 'MarkerSize', 5);
grid on;
xlabel('经度 (Longitude)');
ylabel('纬度 (Latitude)');
title(['Starlink-1008 星下点轨迹 (2圈)']);
axis([-180 180 -90 90]);
hold on;

% 画一条参考线：赤道
yline(0, 'k--');
text(-170, 5, '赤道', 'FontSize', 10);

fprintf('绘图完成！你能看到那个波浪线吗？\n');



%% 接上一步的变量 (请确保你刚才的 workspace 变量还在，或者把这段加在刚才代码的后面)

% ========================================================
% 步骤 5: 设置干扰源 & 计算多普勒 (Doppler)
% ========================================================

% 1. 自动选择干扰源位置 (取轨迹中间的一个点，保证卫星肯定路过)
% 我们取第 200 个采样点作为干扰源位置，这样卫星会在仿真开始后不久飞过它
target_idx = 200; 
lat_target_deg = lat_store(target_idx);
lon_target_deg = lon_store(target_idx);
alt_target = 0; % 假设干扰源在地面，高度为 0

fprintf('干扰源位置已设置: Lat = %.4f°, Lon = %.4f°\n', lat_target_deg, lon_target_deg);

% 2. 坐标系转换准备 (LLA -> ECEF)
% 我们需要把经纬度转成 XYZ 坐标才能算速度和距离
% 定义辅助函数 (简单的球体模型，如果要高精度可以用 wgs84Ellipsoid)
RE = 6371000; % 地球平均半径 (m)

% 干扰源 ECEF 坐标
[x_tgt, y_tgt, z_tgt] = sph2cart(deg2rad(lon_target_deg), deg2rad(lat_target_deg), RE);
pos_target = [x_tgt; y_tgt; z_tgt];

% 3. 计算卫星的 ECEF 坐标序列和速度矢量
% 注意：刚才我们只存了经纬度，现在需要重新算一下 XYZ
pos_sat = zeros(3, length(t));
vel_sat = zeros(3, length(t));

% 为了算速度，我们需要对位置求导，或者利用轨道公式
% 这里用差分法简单估算速度 (v = (p2 - p1) / dt)
% 先算出所有时刻的位置
for k = 1:length(t)
    % 简单的球坐标转直角坐标 (忽略扁率，简化模型)
    % 注意：r 是刚才代码里算的轨道半径 (a * ...) 那个公式里的 r
    % 这里为了简化，我们重新近似算一下
    
    % 重复一下刚才的核心计算以便获取 r 和 lambda, phi
    dt_k = t(k);
    M = M0 + n_rad_s * dt_k;
    E = M; for i=1:10, E=E-(E-ecc*sin(E)-M)/(1-ecc*cos(E)); end
    v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
    r_k = a * (1-ecc^2) / (1+ecc*cos(v)); % 轨道半径
    u = omega + v;
    
    % 惯性系坐标
    x_node = r_k * cos(u);
    y_node = r_k * sin(u);
    
    % 旋转到 ECEF (考虑地球自转 We * dt)
    % 这是一个旋转矩阵变换
    Om = RAAN - We * dt_k; % 修正后的升交点经度
    i = inc;
    
    % 3D 旋转公式 (从轨道面 -> ECEF)
    X = x_node * cos(Om) - y_node * cos(i) * sin(Om);
    Y = x_node * sin(Om) + y_node * cos(i) * cos(Om);
    Z = y_node * sin(i);
    
    pos_sat(:, k) = [X; Y; Z];
end

% 计算速度 (差分法)
for k = 2:length(t)-1
    vel_sat(:, k) = (pos_sat(:, k+1) - pos_sat(:, k-1)) / (2 * (t(2)-t(1)));
end
% 补齐首尾
vel_sat(:, 1) = vel_sat(:, 2);
vel_sat(:, end) = vel_sat(:, end-1);

% 4. 计算多普勒频移 (Doppler Shift)
f_c = 1.57542e9; % 载波频率 (例如 GPS L1: 1.57542 GHz)
c = 299792458;   % 光速
doppler_freq = zeros(1, length(t));
range_dist = zeros(1, length(t));

for k = 1:length(t)
    % 相对位置矢量 (卫星 - 目标)
    r_vec = pos_sat(:, k) - pos_target;
    dist = norm(r_vec);
    range_dist(k) = dist;
    
    % 单位方向矢量 (Line of Sight)
    u_vec = r_vec / dist;
    
    % 相对速度 (假设目标不动，只考虑卫星速度)
    v_vec = vel_sat(:, k);
    
    % 径向速度 (投影)
    v_radial = dot(v_vec, u_vec);
    
    % 多普勒公式: fd = - (v_radial / c) * fc
    doppler_freq(k) = - (v_radial / c) * f_c;
end

% ========================================================
% 6. 绘图：多普勒曲线 (S-Curve)
% ========================================================
figure('Color','w');
subplot(2,1,1);
plot(t, doppler_freq/1000, 'LineWidth', 2);
grid on;
title('多普勒频移曲线 (Doppler S-Curve)');
ylabel('频移 (kHz)');
xlabel('时间 (s)');

subplot(2,1,2);
plot(t, range_dist/1000, 'r', 'LineWidth', 2);
grid on;
title('卫星-干扰源 距离变化');
ylabel('距离 (km)');
xlabel('时间 (s)');

fprintf('多普勒计算完成！\n');

%% 接上一步的代码...

% ========================================================
% 步骤 7: 核心算法 - 单星无源定位 (加噪声 + 最小二乘求解)
% ========================================================

fprintf('\n开始进行定位解算...\n');

% 1. 模拟真实测量 (加入高斯白噪声)
sigma_f = 10; % 测频误差 10Hz (典型的接收机性能)
noise = sigma_f * randn(size(doppler_freq));
f_measured = doppler_freq + noise;

% 选取一段有效数据进行定位 (比如 1500s 到 2500s，卫星过顶这15分钟)
% 实际工程中，我们只在卫星能收到信号时才有数据
valid_idx = 1500:2500; 
t_obs = t(valid_idx);
f_obs = f_measured(valid_idx);
pos_sat_obs = pos_sat(:, valid_idx);
vel_sat_obs = vel_sat(:, valid_idx);

% 2. 定义最小二乘的目标函数
% 我们要找一个位置 x (lat, lon, alt)，使得理论多普勒和测量多普勒的误差平方和最小
% 状态变量 x = [x_ecef; y_ecef; z_ecef]

% 初始猜测 (Initial Guess) - 假设我们只知道目标在地球表面，随便猜一个
% 如果猜得太远，可能不收敛。这里我们故意猜偏一点 (比如真值 + 500km)
x0 = pos_target + [500000; 500000; 500000]; 

% 使用 MATLAB 的 lsqnonlin (非线性最小二乘求解器)
% 匿名函数定义误差: error = f_model - f_meas
cost_func = @(x) doppler_equation(x, pos_sat_obs, vel_sat_obs, f_c, c) - f_obs;

options = optimoptions('lsqnonlin', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt');
[x_est, resnorm] = lsqnonlin(cost_func, x0, [], [], options);

% 3. 结果评估
err_pos = norm(x_est - pos_target);
fprintf('------------------------------------------------\n');
fprintf('真实位置 (ECEF): [%.1f, %.1f, %.1f]\n', pos_target(1), pos_target(2), pos_target(3));
fprintf('解算位置 (ECEF): [%.1f, %.1f, %.1f]\n', x_est(1), x_est(2), x_est(3));
fprintf('定位误差: %.2f 米\n', err_pos);
fprintf('------------------------------------------------\n');

% --------------------------------------------------------
% 辅助函数: 多普勒方程
% --------------------------------------------------------
function f_d = doppler_equation(x_target, pos_s, vel_s, f_c, c)
    % x_target: 3x1 猜测的目标位置
    % pos_s: 3xN 卫星位置
    % vel_s: 3xN 卫星速度
    
    N = size(pos_s, 2);
    f_d = zeros(1, N);
    
    for k = 1:N
        r_vec = pos_s(:, k) - x_target;
        dist = norm(r_vec);
        u_vec = r_vec / dist;
        v_vec = vel_s(:, k);
        v_radial = dot(v_vec, u_vec);
        f_d(k) = - (v_radial / c) * f_c;
    end
end