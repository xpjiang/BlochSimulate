%% -------- 单点 FID 仿真 --------
clear; clc;

%% -------- 手动生成 2 个自旋 --------
% 磁化矢量 [3 x N]
Spins.M = [ ...
    0   0;   % Mx
    0   0;   % My
    1   1];  % Mz
% 空间位置 [3 x N]（单位：m）
Spins.Pos = [ ...
     -1e-3      1e-3;   % x
     0      0;      % y
     0      0];     % z
% 弛豫参数 [1 x N]
Spins.T1 = [1.5  1.5];   % s
Spins.T2 = [0.08 0.08];  % s
% 质子密度（M0）
Spins.PD = [1.0  0.8];   % 第二个自旋信号弱一点
% 静磁场不均匀（Hz）
Spins.dB = [0  0];      % 第二个自旋 off-resonance 50 Hz


dt = 1e-6;                   % 1 us
total_time = 300e-3;
MTime = round(total_time / dt);
gamma_rad = 2*pi*42.57e6;

% --- RF事件 ---
tPulse = 5e-3; %100us
tPulseStart = 1e-3; %100us
tPulseEnd = tPulseStart + tPulse;%s

flip = pi/2;                  % 90°
B1_amp = flip / (gamma_rad * tPulse);

RF.B1 = zeros(1, MTime);
RF.df = zeros(1, MTime);
RF.is_active = false(1, MTime);

RF.B1(round(tPulseStart/dt):round(tPulseEnd/dt)) = B1_amp;  % 第 1 个时间步打 RF
RF.is_active(round(tPulseStart/dt):round(tPulseEnd/dt)) = true;

%  --- 梯度事件 ---
tGrad = 5e-3; %100us
tGradStart = tPulseEnd + 100e-3; %100us
tGradEnd = total_time;%s

Grads.Gxyz = zeros(3, MTime);
Grads.Gxyz(1,round(tGradStart/dt):round(tGradEnd/dt)) = 0.8e-3; % 0.08 mT/m
Grads.is_active = false(1, MTime);
Grads.is_active(round(tGradStart/dt):round(tGradEnd/dt)) = true;

%Spins.dB  = -0.08e-3 * 0.01 * 42.57e6; % Hz（off-resonance）

[M_final, M_hist] = bloch_evolve_active_indexed( ...
    Spins, RF, Grads, dt, true);

full_signal = squeeze(sum(M_hist(1,:,:) + 1i*M_hist(2,:,:), 2));

t_ms = (1:MTime) * dt * 1000;
plot_mri_sequence(Spins, RF, Grads, t_ms, full_signal);

