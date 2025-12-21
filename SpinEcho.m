%% 球形分布下的 Spin Echo 仿真
clear; clc;

% --- 1. 物理环境与自旋位置初始化 ---
N_spins = 1000;
R_ball = 0.05; % 50mm 半径 = 100mm 直径

% 生成球内均匀分布的坐标
r = R_ball * (rand(1, N_spins).^(1/3));      % 半径补偿
phi = 2 * pi * rand(1, N_spins);             % 方位角
theta = acos(2 * rand(1, N_spins) - 1);      % 极角

% 转换为笛卡尔坐标 [3 x N]
Pos = [r .* sin(theta) .* cos(phi); ...
       r .* sin(theta) .* sin(phi); ...
       r .* cos(theta)];

% --- 2. 物理参数与场不均匀性 ---
% 模拟一个线性磁场梯度偏差 dB (模拟真实环境中的不完美匀场)
% 假设沿 z 轴存在 10Hz/mm 的场偏差
dB_inhomogeneity = 0 * r; % Hz
T2 = 0.08 * ones(1, N_spins);
PD = ones(1, N_spins);
M_init = [zeros(2, N_spins); PD];

% --- 3. 序列设置 (1us 步长) ---
dt = 1e-6;
TE = 20e-3;
total_time = 60e-3;
MTime = round(total_time / dt);


B1_seq = zeros(1, MTime);
Gx_seq = zeros(1, MTime);
Gx_seq(:) = 0.00e-3; % 0.08 mT/m
is_sampling = false(1, MTime);

% 定义脉冲与采样 (复用之前的逻辑)
gamma_rad = 42.57e6 * 2 * pi;
T90 = 1e-3; %1ms
T180 = 1e-3; %1ms
B1_seq(5:round(T90/dt + 5)) = (pi/4) / (gamma_rad * T90); % 90 deg
idx_180 = round((TE/2-T180/2)/dt) : round((TE/2+T180/2)/dt);
B1_seq(idx_180) = 1i * (3*pi/4) / (gamma_rad * T180); % 180 deg

Spins.M = M_init;
Spins.Pos = Pos;
Spins.T1 = 1.0*ones(1,N_spins);
Spins.T2 = T2;
Spins.PD = PD;
Spins.dB = dB_inhomogeneity;

RF.B1 = B1_seq;
RF.df = zeros(1, MTime);
RF.df = any(abs(RF.B1) > 1e-9, 1) * 20000;
RF.is_active = abs(RF.B1) > 1e-12;

Grads.Gxyz = [Gx_seq; Gx_seq;Gx_seq];
Grads.is_active = any(abs(Grads.Gxyz) > 1e-9, 1);

% --- 4. 运行演化 ---
[~, M_history] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);

full_signal = squeeze(sum(M_history(1,:,:) + 1i*M_history(2,:,:), 2)) / N_spins;

t_ms = (1:MTime) * dt * 1000;
plot_mri_sequence(Spins, RF, Grads, t_ms, full_signal);