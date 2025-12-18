%% 从 Mz 开始的完整布洛赫仿真 (激发 + 采样)
clear; clc; close all;

% --- 1. 基础物理参数 ---
gamma = 42.57e6;            % Hz/T
gamma_rad = gamma * 2 * pi;
T1 = 0.20;                   % 1s
T2 = 0.01;                  % 50ms
M0 = 1;

% --- 2. 模拟环境设置 (考虑频率偏移) ---
df_offset = 100;            % 100Hz 的频率偏移 (Off-resonance)
dt = 1e-5;                  % 时间步长 (10us)，激发过程需要更细的时间步长

% --- 3. 阶段一：RF 激发 (RF Excitation) ---
flip_angle = 90;            % 翻转角 (度)
t_rf = 0.04;                % RF 脉冲持续时间 1ms
% 计算达到该翻转角所需的 B1 强度: flip_angle = gamma_rad * B1 * t_rf
B1_amp = (deg2rad(flip_angle)) / (gamma_rad * t_rf); 

time_rf = 0:dt:t_rf;
M = [0; 0; M0];             % 初始状态：处于 Mz 轴
M_exc_history = zeros(3, length(time_rf));

for t = 1:length(time_rf)
% A. 构建 3x3 算子矩阵 A
    % 旋转部分 (假设 RF 沿 x 轴，考虑 10Hz 的微小失谐)
    df = 10; 
    omega_x = gamma_rad * B1_amp;
    omega_z = 2 * pi * df;
    
    % 弛豫部分
    R1 = 1/T1;
    R2 = 1/T2;
    
    A = [ -R2,  -omega_z,   0;
           omega_z,  -R2,   -omega_x;
           0,        omega_x,  -R1 ];
       
    % B. 构建扩展矩阵 F (4x4)
    % dM/dt = A*M + [0;0;M0/T1] -> d[M;1]/dt = [A, C; 0, 0] * [M;1]
    C = [0; 0; M0/T1];
    F = zeros(4,4);
    F(1:3, 1:3) = A;
    F(1:3, 4) = C;
    
    % C. 矩阵指数求解 (离散步长演化)
    % 在 MATLAB 中，expm 是求解矩阵指数的标准方法
    ExpFdt = expm(F * dt);
    
    % D. 更新状态
    M_aug = ExpFdt * [M; 1];
    M = M_aug(1:3);
    M_exc_history(:, t) = M;
end

% --- 4. 阶段二：自由感应衰减 (FID Sampling) ---
t_fid = 0.1;                % 采样 100ms
time_fid = 0:dt:t_fid;
M_fid_history = zeros(3, length(time_fid));

% 此时 B1 = 0，仅剩频率偏移导致的进动和 T1/T2 弛豫
for t = 1:length(time_fid)
    % 进动角速度
    omega_z = 2 * pi * df_offset;
    
    % 构造旋转 + 弛豫核
    phi = omega_z * dt;
    E1 = exp(-dt/T1);
    E2 = exp(-dt/T2);
    
    % 矩阵更新
    R_z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    M = R_z * M;                     % 进动
    M = [E2 0 0; 0 E2 0; 0 0 E1] * M + [0; 0; M0*(1-E1)]; % 弛豫
    
    M_fid_history(:, t) = M;
end

% --- 5. 结果可视化 ---
figure('Color', 'w', 'Position', [100 100 900 400]);

% 绘制磁化矢量的轨迹
subplot(1,2,1);
t_total = [time_rf, time_rf(end) + time_fid];
M_total = [M_exc_history, M_fid_history];
plot(t_total*1000, M_total(1,:), 'r', 'DisplayName', 'Mx'); hold on;
plot(t_total*1000, M_total(2,:), 'g', 'DisplayName', 'My');
plot(t_total*1000, sqrt(M_total(2,:).*M_total(2,:) + M_total(1,:).*M_total(1,:)),'DisplayName', 'Mxy');
plot(t_total*1000, M_total(3,:), 'b', 'DisplayName', 'Mz', 'LineWidth', 2);
xline(t_rf*1000, '--', 'RF End');
xlabel('Time (ms)'); ylabel('Magnetization');
title('Magnetization Evolution (Excitation + FID)');
legend; grid on;

% 绘制横截面 Mxy 的演变
subplot(1,2,2);
plot(M_total(1,:), M_total(2,:));
xlabel('Mx'); ylabel('My');
title('Mxy Plane Trajectory');
axis equal; grid on;