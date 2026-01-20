%% MRI 选层激发仿真：ADC 采样与 Slice Profile 重建
clear; clc; close all;

% --- 1. 协议参数与物理常量 ---
slice_thickness = 0.01;      % 10mm
FOV_z = 0.10;                % 100mm
BW_adc = 20000;              % 20kHz
N_samples = 128;             
gamma = 42.57e6;             
gamma_rad = 2 * pi * gamma;
dt = 1e-6;                   % 仿真步长
t_ramp = 0.0005; 
t_flat_ss = 0.003;

% --- 2. 构造 Spins 结构体 (R2 分布) ---
nVoxel = 256;               % 增加体素密度提高重绘质量
nSpinsPerVoxel = 24;        
voxel_size = FOV_z / nVoxel;
z_voxel_centers = linspace(-FOV_z/2, FOV_z/2, nVoxel);
phi_golden = (1 + sqrt(5)) / 2;
alpha_r2 = 1 / phi_golden;
r2_seq = mod(0.5 + (0:nSpinsPerVoxel-1) * alpha_r2, 1) - 0.5;
z_micro_offsets = r2_seq * voxel_size;
Z_all = reshape(z_voxel_centers' + z_micro_offsets, 1, []);
num_spins = length(Z_all);

Spins.M   = repmat([0; 0; 1], 1, num_spins);
Spins.Pos = [zeros(2, num_spins); Z_all]; 
Spins.T1  = ones(1, num_spins) * 1.0;     
Spins.T2  = ones(1, num_spins) * 0.1;     
Spins.PD  = ones(1, num_spins);
Spins.dB  = (rand(1, num_spins) - 0.5) * 5; % 略微减小 dB 观察纯净 Profile

% --- 3. 构造 RF 和 Grads 结构体 ---
% A. 选层阶段
G_ss_amp = G_ss_amp_calc();
G_ss_wave = [linspace(0, G_ss_amp, t_ramp/dt), ones(1, t_flat_ss/dt)*G_ss_amp, linspace(G_ss_amp, 0, t_ramp/dt)];

% B1 脉冲
rf_wave = zeros(size(G_ss_wave));
t_ss_axis = (0:length(G_ss_wave)-1)*dt;
rf_mask = (t_ss_axis >= t_ramp) & (t_ss_axis <= t_ramp + t_flat_ss);
t_rf_rel = t_ss_axis(rf_mask) - (t_ramp + t_flat_ss/2);
raw_sinc = sinc((4/t_flat_ss) * t_rf_rel);
B1_peak = deg2rad(90) / (gamma_rad * sum(raw_sinc)*dt);
rf_wave(rf_mask) = B1_peak * raw_sinc;

% B. 重聚与预相位 (Rewinder)
t_flat_pre = 0.002;
G_pre_amp = calculate_pre_amp();
G_pre_wave = [linspace(0, G_pre_amp, t_ramp/dt), ones(1, ceil(t_flat_pre/dt))*G_pre_amp, linspace(G_pre_amp, 0, t_ramp/dt)];

% C. 读取阶段 (ADC 采样在此阶段平顶发生)
t_read_flat = N_samples / BW_adc;
G_read_amp = BW_adc / (gamma * FOV_z);
G_read_wave = [linspace(0, G_read_amp, t_ramp/dt), ones(1, ceil(t_read_flat/dt))*G_read_amp, linspace(G_read_amp, 0, t_ramp/dt)];

% 拼接总波形
Gz_total = [G_ss_wave, G_pre_wave, G_read_wave];
B1_total = [rf_wave, zeros(1, length(G_pre_wave)), zeros(1, length(G_read_wave))];

% 确定 ADC 采样窗索引
adc_start_idx = length(G_ss_wave) + length(G_pre_wave) + round(t_ramp/dt);
adc_end_idx = adc_start_idx + round(t_read_flat/dt) - 1;

RF.B1 = B1_total;
RF.df = zeros(size(B1_total));
RF.is_active = (abs(B1_total) > 0);
RF.is_adc = zeros(size(B1_total));
RF.is_adc(adc_start_idx:adc_end_idx) = 1; % 标记 ADC 窗口

Grads.Gxyz = [zeros(2, length(Gz_total)); Gz_total];
Grads.is_active = (abs(Gz_total) > 0);

% --- 4. 执行演化 ---
[~, M_history] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);
signal = squeeze(sum(M_history(1,:,:) + 1i*M_history(2,:,:), 2)) / num_spins;

% --- 5. ADC 采样过程 ---
t_total = (0:length(signal)-1) * dt;
% 理想采样时间点 (居中采样)
t_adc_samples = (0:N_samples-1) / BW_adc + t_total(adc_start_idx);
% 从仿真信号中插值得到 ADC 原始数据
adc_data = interp1(t_total, signal, t_adc_samples, 'linear', 'extrap');

% --- 6. IFT 重建层剖面 ---
recon_profile = fftshift(ifft(ifftshift(adc_data)));
% 计算重构坐标轴: z = f / (gamma * G_read)
% 频率分辨率 df = BW_adc / N_samples
dz = BW_adc / (gamma * G_read_amp * N_samples);
z_axis_recon = ((0:N_samples-1) - N_samples/2) * dz;

% --- 7. 绘图与验证 ---
t_ms = t_total * 1000;
% 调用你改进后的 plot_mri_sequence (包含三轴梯度和 ADC 背景)
ax_handles = plot_mri_sequence2(Spins, RF, Grads, t_ms, signal);

% 创建专门的验证窗口
figure('Color', 'w', 'Name', 'ADC Signal & Slice Profile Verification');

% 子图 A: ADC 采样到的时域信号 (k-space)
subplot(2,1,1);
t_adc_ms = (t_adc_samples - t_adc_samples(1)) * 1000; % 采样相对时间
plot(t_adc_ms, abs(adc_data), 'k', 'LineWidth', 2, 'DisplayName', 'Magnitude'); hold on;
plot(t_adc_ms, real(adc_data), 'b', 'DisplayName', 'Real (I)');
plot(t_adc_ms, imag(adc_data), 'r', 'DisplayName', 'Imag (Q)');
grid on; xlabel('ADC Sampling Time (ms)'); ylabel('Signal Amplitude');
title('Figure 2a: Raw ADC Sampled Signal (Time Domain)');
legend('Location', 'northeast');

% 子图 B: 重建后的层剖面 (Image Domain)
subplot(2,1,2);
% 归一化处理方便对比
norm_profile = abs(recon_profile) / max(abs(recon_profile));
plot(z_axis_recon * 1000, norm_profile, 'LineWidth', 2, 'Color', [0 0.5 0]);
hold on;
% 绘制理论层厚边界
y_lims = [0 1.2];
patch([-slice_thickness/2, slice_thickness/2, slice_thickness/2, -slice_thickness/2]*1000, ...
      [0 0 1.1 1.1], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'r', 'LineStyle', '--', 'DisplayName', 'Target Slice');
grid on; 
xlabel('Position z (mm)'); ylabel('Normalized Intensity');
title('Figure 2b: Reconstructed Slice Profile (IFT of ADC Data)');
legend('Reconstructed', 'Theoretical Target');
xlim([-FOV_z/2, FOV_z/2] * 1000); % 限制在 FOV 范围内观察
ylim(y_lims);

%% ===== 辅助函数 =====
function val = G_ss_amp_calc()
    val = (4/0.003) / (42.57e6 * 0.01);
end
function val = calculate_pre_amp()
    G_ss_a = G_ss_amp_calc();
    G_rd_a = 20000 / (42.57e6 * 0.10);
    area_rew = G_ss_a * (0.003/2 + 0.0005/2);
    area_pre = G_rd_a * (128/2/20000 + 0.0005/2);
    val = -(area_rew + area_pre) / (0.002 + 0.0005);
end