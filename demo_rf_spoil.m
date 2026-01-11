%% Steady State Approach Simulation (Mxy at TE vs. TR Number)
% 目的: 绘制每个 TR 在 TE 时刻的信号，观察稳态建立过程
% Author: MRI Expert (Gemini)

clear; close all; clc;

% --- 1. 参数设置 ---
nTR = 200;              % 模拟 200 个 TR
TR = 5e-3;             % TR = 5 ms
TE = 1.5e-3;            % TE = 1.5 ms
FlipAngle = 26;         % 翻转角
dt = 10e-6;             % 仿真步长

% --- 2. 自旋系统设置 ---
% 必须引入频率分散 (dB) 来模拟 Gradient Spoil 的效果
nSpins = 200;
Spins.Pos = zeros(3, nSpins); 
Spins.T1  = 200e-3 * ones(1, nSpins);  % T1 = 200ms
Spins.T2  = 80e-3  * ones(1, nSpins);  % T2 = 80ms
Spins.PD  = ones(1, nSpins);
% 关键：模拟体素内的频率分布 (Gradient Dephasing)
Spins.dB  = linspace(-60, 60, nSpins); 

% 初始磁化状态 (M0)
M_init = [zeros(2, nSpins); ones(1, nSpins)];

% 预计算参数
n_steps_TR = round(TR/dt);
idx_TE     = round(TE/dt); % TE 在数组中的索引
gamma = 42.57e6;
t_rf = 1e-3; % 1ms
B1_amp = (FlipAngle/360) / (gamma * t_rf);

%% --- 仿真 1: 无 RF 扰相 (Standard GRE) ---
fprintf('Simulating Standard GRE (No RF Spoil)...\n');
Spins.M = M_init;
Sig_NoSpoil = zeros(1, nTR);

for n = 1:nTR
    % 相位恒定为 0
    RF_phase = 0;
    
    % 生成序列 & 演化
    [RF, Grads] = make_sequence(n_steps_TR, B1_amp, RF_phase, dt);
    [M_final, M_hist] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);
    
    % 提取 TE 时刻信号 (矢量和)
    M_at_TE = sum(M_hist(:,:,idx_TE), 2);
    Sig_NoSpoil(n) = M_at_TE(1) + 1i*M_at_TE(2);
    
    % 更新状态
    Spins.M = M_final;
end

%% --- 仿真 2: 有 RF 扰相 (RF Spoiled / SPGR) ---
fprintf('Simulating RF Spoiled GRE (SPGR)...\n');
Spins.M = M_init;
Sig_RFSpoil = zeros(1, nTR);

deg2rad = pi/180;

rf_inc=0;       % 累积相位
rfSpoilingInc = 117;      % 增量步长 (117度)
rf_phase=0;
for n = 1:nTR
    % 二次相位更新: Phi(n) = Phi(n-1) + n * 117
    phaseOffset = rf_phase * deg2rad;
    
    % 生成序列 & 演化
    [RF, Grads] = make_sequence(n_steps_TR, B1_amp, phaseOffset, dt);
    [M_final, M_hist] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, true);
    
    % 提取 TE 时刻信号
    M_at_TE = sum(M_hist(:,:,idx_TE), 2);
    % 解调 (Demodulation): 移除发射相位的影响
    raw_sig = M_at_TE(1) + 1i*M_at_TE(2);
    Sig_RFSpoil(n) = raw_sig * exp(-1i * phaseOffset);
    
    % 更新相位
    % rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    % rf_phase=mod(rf_phase+rf_inc, 360.0);
    rf_phase= mod(rf_phase + rfSpoilingInc, 360.0); %% why
    Spins.M = M_final;
end

%% --- 3. 结果绘制 ---
figure('Color','w', 'Position', [200 200 900 600]);

% 子图 1: 幅度演变 (Magnitude)
subplot(2,1,1);
plot(1:nTR, abs(Sig_NoSpoil), 'r.-', 'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(1:nTR, abs(Sig_RFSpoil), 'b.-', 'LineWidth', 1, 'MarkerSize', 8);
grid on;
title('Approach to Steady State (Magnitude at TE)');
xlabel('TR Number'); ylabel('|M_{xy}| Signal');
legend('No RF Spoil (Oscillatory/SSFP)', 'RF Spoiled (Smooth/T1-weighted)', 'Location', 'best');
% 添加理论稳态值的注释
text(nTR/2, mean(abs(Sig_NoSpoil(end-20:end))), '\leftarrow Higher Signal (T2/T1 dependent)', 'Color','r');
text(nTR/2, mean(abs(Sig_RFSpoil(end-20:end))), '\leftarrow Pure T1 Signal (Ernst)', 'Color','b','VerticalAlignment','top');

% 子图 2: 相位演变 (Phase)
subplot(2,1,2);
plot(1:nTR, angle(Sig_NoSpoil), 'r.'); hold on;
plot(1:nTR, angle(Sig_RFSpoil), 'b.');
grid on;
title('Signal Phase at TE (Demodulated)');
xlabel('TR Number'); ylabel('Phase (rad)');
ylim([-pi pi]);
legend('No RF Spoil', 'RF Spoiled');

sgtitle(['Steady State Analysis: TR=' num2str(TR*1e3) 'ms, \alpha=' num2str(FlipAngle) '^\circ']);

%% ===== 辅助函数 1: 序列生成 =====
function [RF, Grads] = make_sequence(n_steps, B1_val, phase, dt)
    RF.B1 = zeros(1, n_steps);
    RF.df = zeros(1, n_steps);
    RF.is_active = false(1, n_steps);
    Grads.Gxyz = zeros(3, n_steps);
    Grads.is_active = false(1, n_steps);
    
    % RF 脉冲 (1ms)
    n_rf = round(1e-3/dt);
    RF.B1(1:n_rf) = B1_val * exp(1i * phase);
    RF.is_active(1:n_rf) = true;
    
    % Spoiler Gradient (TR 最后 1ms)
    % spoiler_dur = round(1e-3 / dt);
    % start_spoil = n_steps - spoiler_dur + 1;
    % Grads.Gxyz(3, start_spoil:end) = 0.02; % Z轴梯度
    % Grads.is_active(start_spoil:end) = true;
end

