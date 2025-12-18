%% MRI 选层激发与梯度回波读取完整仿真
clear; clc; close all;

% --- 1. 协议参数 (User Inputs) ---
slice_thickness = 0.01;      % 10mm
FOV_z = 0.10;                 % 100mm
BW_adc = 20000;              % 20kHz ADC带宽
N_samples = 128;             % 采样点数
gamma = 42.57e6;             % Hz/T

t_ramp = 0.0005; t_flat_ss = 0.003;

% --- 2. 物理参数计算 ---
% --- 射频脉冲强度计算 ---
flip_angle_deg = 90;               % 目标翻转角 (度)
flip_angle_rad = deg2rad(flip_angle_deg);
% 计算 Sinc 脉冲的积分因子 (对于标准 sinc，积分面积约为 1/BW_rf)
% 精确计算：B1_peak = flip_angle / (gamma * area_under_sinc)
% 在离散仿真中，我们直接用 sum(B1) * dt = flip_angle / gamma 来归一化

TBP = 4;                             % Time-Bandwidth Product (典型值)
t_rf = t_flat_ss;                    % 脉冲持续时间 (3ms)
BW_rf = TBP / t_rf;                  % 自动计算带宽: 4 / 0.003 = 1333.3 Hz
% --- 对应的梯度强度计算 ---
G_ss_amp = BW_rf / (gamma * slice_thickness);
G_read_amp = BW_adc / (gamma * FOV_z);         % 读取梯度强度 (~4.7 mT/m)
dt = 1e-5;                   % 仿真步长 (10us)



% --- 3. 构造时序波形 ---
% (A) 选层阶段 (Slice Selection)

G_ss = [linspace(0, G_ss_amp, t_ramp/dt), ones(1, t_flat_ss/dt)*G_ss_amp, linspace(G_ss_amp, 0, t_ramp/dt)];

% (B) 重聚与预相位 (Rewinder + Prephase)
% 面积 = 选层补偿面积 + 读取预相位面积 (为了让回波在采样窗中心)
t_flat_pre = 0.0015;
area_rew = G_ss_amp * (t_flat_ss/2 + t_ramp/2);
area_pre = G_read_amp * (N_samples/2 / BW_adc + t_ramp/2);
G_pre_amp = -(area_rew + area_pre) / (t_flat_pre + t_ramp);
G_pre = [linspace(0, G_pre_amp, t_ramp/dt), ones(1, t_flat_pre/dt)*G_pre_amp, linspace(G_pre_amp, 0, t_ramp/dt)];

% (C) 读取阶段 (Readout)
t_read_flat = N_samples / BW_adc; % 6.4ms
G_read = [linspace(0, G_read_amp, t_ramp/dt), ones(1, t_read_flat/dt)*G_read_amp, linspace(G_read_amp, 0, t_ramp/dt)];

% 拼接总梯度
Gz_total = [G_ss, G_pre, G_read];
t_total = (0:length(Gz_total)-1) * dt;

% 构造RF脉冲 (中心对准SS平顶中心)
B1 = zeros(size(t_total));
t_center_ss = t_ramp + t_flat_ss/2;
rf_mask = (t_total >= t_ramp) & (t_total <= t_ramp + t_flat_ss);

% 1. 先生成原始 Sinc 波形 (无单位)
raw_sinc = sinc(BW_rf * (t_total(rf_mask) - t_center_ss));

% 2. 计算归一化系数，确保积分等于目标翻转角
% 这里的 gamma 是 Hz/T，所以需要乘以 2*pi 转换成 rad/(s·T)
gamma_rad = gamma * 2 * pi; 
current_area = sum(raw_sinc) * dt; % 计算当前离散脉冲的面积
B1_peak = flip_angle_rad / (gamma_rad * current_area); % 核心公式

% 3. 赋予物理单位 (Tesla)
B1(rf_mask) = B1_peak * raw_sinc;


% 确定ADC采集区间 (Readout平顶期)
adc_start_time = t_ramp + t_flat_ss + t_ramp + t_ramp + t_flat_pre + t_ramp + t_ramp;
adc_end_time = adc_start_time + t_read_flat;
adc_mask = (t_total >= adc_start_time) & (t_total <= adc_end_time);

% --- 4. 多自旋布洛赫仿真 (增加 B0 不均匀性) ---
N_spins = 512;
z_pos = linspace(-3*FOV_z/2, 3*FOV_z/2, N_spins);
spins = NMRSpin.empty(N_spins, 0);

% 定义不均匀性参数
dB0_gradient = 5; % 线性偏磁场 (Hz/m)，模拟磁场线性漂移
dB0_random = 2;     % 随机偏磁场 (Hz)，模拟微观局部不均匀性

for i = 1:N_spins
    % 计算该位置的局部磁场偏差
    % delta_B = 线性项 + 随机项
    local_dB = (z_pos(i) * dB0_gradient) + (randn() * dB0_random);
    
    % 设置质子密度：10mm 水模
    % pd = (abs(z_pos(i)) <= slice_thickness/2);
    pd = 1;
    
    % 创建自旋对象，传入 local_dB 作为第四个参数 (dB)
    spins(i) = NMRSpin(pd, 1.0, 0.01, local_dB, [0, 0, z_pos(i)]);
end

S_t = zeros(size(t_total));
for t = 1:length(t_total)
    M_sum = 0;
    for i = 1:N_spins
        spins(i).evolve(dt, B1(t), [0, 0, Gz_total(t)], 0);
        M_sum = M_sum + (spins(i).M(1) + 1i*spins(i).M(2));
    end
    S_t(t) = M_sum;
end

% --- 5. ADC 采样与重建 (真正实现 BW_adc 采样) ---
dt_adc = 1 / BW_adc; 
t_adc_samples = adc_start_time : dt_adc : (adc_start_time + t_read_flat - dt_adc);
adc_signal_real_data = interp1(t_total, S_t, t_adc_samples, 'linear');

% 计算重建
recon = fftshift(fft(fftshift(adc_signal_real_data)));
dist_axis = linspace(-FOV_z/2, FOV_z/2, length(recon)) * 1000;

% --- 6. 综合绘图展示 ---
figure('Color', 'w', 'Position', [100 50 1100 900]);

% A. 序列时序与 ADC 窗口
subplot(4,1,1);
hold on; yyaxis left;
plot(t_total*1000, Gz_total*1000, 'b', 'LineWidth', 2);
patch_t = [adc_start_time, adc_end_time, adc_end_time, adc_start_time]*1000;
patch_g = [min(Gz_total*1000), min(Gz_total*1000), max(Gz_total*1000), max(Gz_total*1000)];
fill(patch_t, patch_g, 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
ylabel('Gz (mT/m)'); yyaxis right;
plot(t_total*1000, real(B1)*1e6, 'r'); ylabel('B1 (\muT)');
title('1. Sequence Layout (ADC Window in Yellow)'); grid on;

% B. 连续信号演化 (整个过程)
subplot(4,1,2);
plot(t_total*1000, abs(S_t), 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
plot(t_adc_samples*1000, abs(adc_signal_real_data), 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
title('2. Transverse Magnetization Evolution & ADC Sampling Points');
ylabel('Signal'); legend('Continuous', 'ADC Samples'); grid on;

% C. ADC 采样时域复信号 (Raw k-space Data)
subplot(4,1,3);
t_relative = (0:length(adc_signal_real_data)-1) * dt_adc * 1000; % 从0开始的采样时间
plot(t_relative, real(adc_signal_real_data), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Real (I)'); hold on;
plot(t_relative, imag(adc_signal_real_data), 'r-', 'LineWidth', 1.2, 'DisplayName', 'Imag (Q)');
plot(t_relative, abs(adc_signal_real_data), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Magnitude');
title('3. ADC Raw Time-Domain Signal (Detected Echo)');
xlabel('Sampling Time (ms)'); ylabel('Amplitude'); legend; grid on;

% D. 1D 重建剖面
subplot(4,1,4);
plot(dist_axis, abs(recon), 'LineWidth', 2, 'Color', [0 0.5 0]);
title('4. Final Reconstructed 1D Profile');
xlabel('Z Position (mm)'); ylabel('Intensity');
xline(-slice_thickness/2*1000, '--r'); xline(slice_thickness/2*1000, '--r');
grid on; xlim([-FOV_z*0.6, FOV_z*0.6]*1000);