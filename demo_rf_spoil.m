%% Steady State Approach Simulation: R2 Distribution & Gradient Spoiling
% Author: MRI Expert (Gemini)
%clear; close all; clc;

% --- 1. 参数设置 ---
nTR = 200;              % 模拟 200 个 TR
TR = 5e-3;              % 5 ms
TE = 1.5e-3;            % 1.5 ms
FlipAngle = 26;         % 翻转角
dt = 10e-6;             % 仿真步长

% --- 2. 自旋系统设置 (R2 Low-discrepancy) ---
nSpins = 2000;
phi_golden = (1 + sqrt(5)) / 2;
alpha_r2 = 1 / phi_golden;
seed = 0.5; 
% 在 Z 轴方向生成 [ -0.5, 0.5 ] 之间的低偏差分布位置
r2_pos = mod(seed + (0:nSpins-1) * alpha_r2, 1) - 0.5;

Spins.Pos = zeros(3, nSpins); 
Spins.Pos(3,:) = r2_pos; 
Spins.T1  = 200e-3 * ones(1, nSpins);  
Spins.T2  = 80e-3  * ones(1, nSpins);  
Spins.PD  = ones(1, nSpins);
Spins.dB  = zeros(1, nSpins);

M_init = [zeros(2, nSpins); ones(1, nSpins)];

% 预计算参数
n_steps_TR = round(TR/dt);
idx_TE     = round(TE/dt); 
gamma = 42.57e6; 
gamma_rad = gamma * 2 * pi;
t_rf = 1e-3; 
B1_amp = (FlipAngle/360) / (gamma * t_rf);

% 梯度强度设计：在 1ms 内产生 2*pi 的体素内散相 (假设体素宽度 1mm)
voxel_width = 1e-3; 
G_amp = 1 / (gamma * 1e-3 * voxel_width); 

%% --- 仿真 1: 无 RF 扰相 (Standard GRE / SSFP) ---
fprintf('Simulating Standard GRE...\n');
Spins.M = M_init;
Sig_NoSpoil = zeros(1, nTR);
for n = 1:nTR
    [RF, Grads] = make_sequence(n_steps_TR, B1_amp, 0, dt, G_amp);
    [M_final, M_hist] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);
    M_at_TE = sum(M_hist(:,:,idx_TE), 2);
    Sig_NoSpoil(n) = (M_at_TE(1) + 1i*M_at_TE(2)) / nSpins;
    Spins.M = M_final;
end

%% --- 仿真 2: 有 RF 扰相 (RF Spoiled / SPGR) ---
fprintf('Simulating RF Spoiled GRE...\n');
Spins.M = M_init;
Sig_RFSpoil = zeros(1, nTR);
deg2rad = pi/180;
rf_inc = 0; rf_phase = 0; rfSpoilingInc = 117;

for n = 1:nTR
    phase_rad = rf_phase * deg2rad;
    [RF, Grads] = make_sequence(n_steps_TR, B1_amp, phase_rad, dt, G_amp);
    [M_final, M_hist] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);
    
    M_at_TE = sum(M_hist(:,:,idx_TE), 2);
    raw_sig = (M_at_TE(1) + 1i*M_at_TE(2)) / nSpins;
    % 解调
    Sig_RFSpoil(n) = raw_sig * exp(-1i * phase_rad);
    
    % 更新相位
    rf_inc = mod(rf_inc + rfSpoilingInc, 360.0);
    rf_phase = mod(rf_phase + rf_inc, 360.0);
    Spins.M = M_final;
end

%% --- 3. 理论分析与绘图 ---
flip_rad = FlipAngle * pi/180;
T1_val = Spins.T1(1); T2_val = Spins.T2(1);
% Ernst Angle 理论公式
E1 = exp(-TR/T1_val);
S_theoretical = (sin(flip_rad) * (1-E1) / (1 - cos(flip_rad)*E1)) * exp(-TE/T2_val);

figure('Color','w', 'Position', [100 100 1000 700]);

% 子图 1: 信号幅度对比
subplot(2,1,1);
plot(1:nTR, abs(Sig_NoSpoil), 'r-', 'LineWidth', 1.2); hold on;
plot(1:nTR, abs(Sig_RFSpoil), 'b-', 'LineWidth', 1.8);
yline(S_theoretical, 'g--', 'Theory (Ernst)', 'LineWidth', 2);
grid on;
title(['Signal Magnitude: Approach to Steady State (\phi_{inc} = ', num2str(rfSpoilingInc), '^{\circ})']);
xlabel('TR Number'); ylabel('|M_{xy}| Signal');
legend('Standard GRE (No RF Spoil)', 'RF Spoiled (SPGR)', 'Ernst Limit');

% 子图 2: 解调后的相位对比
subplot(2,1,2);
plot(1:nTR, angle(Sig_NoSpoil), 'r.'); hold on;
plot(1:nTR, angle(Sig_RFSpoil), 'b.');
yline(0, 'k--');
grid on; grid minor;
title('Demodulated Signal Phase at TE');
xlabel('TR Number'); ylabel('Phase (rad)');
ylim([-pi pi]);
legend('No RF Spoil', 'RF Spoiled');

% 打印误差分析
sim_val = mean(abs(Sig_RFSpoil(end-10:end)));
fprintf('\n--- Analysis ---\n');
fprintf('Ernst Theory: %.6f\n', S_theoretical);
fprintf('RF Spoil Sim: %.6f\n', sim_val);
fprintf('Error: %.2f%%\n', abs(sim_val - S_theoretical)/S_theoretical*100);

%% ===== 辅助函数集 =====
function [RF, Grads] = make_sequence(n_steps, B1_val, phase, dt, G_val)
    RF.B1 = zeros(1, n_steps); RF.df = zeros(1, n_steps); RF.is_active = false(1, n_steps);
    Grads.Gxyz = zeros(3, n_steps); Grads.is_active = false(1, n_steps);
    n_rf = round(1e-3/dt);
    RF.B1(1:n_rf) = B1_val * exp(1i * phase);
    RF.is_active(1:n_rf) = true;
    % Spoiler Gradient
    spoiler_dur = round(1e-3 / dt);
    start_spoil = n_steps - spoiler_dur + 1;
    Grads.Gxyz(3, start_spoil:end) = G_val;
    Grads.is_active(start_spoil:end) = true;
end

