function [M_final, M_history] = bloch_evolve_active_indexed(Spins, RF, Grads, dt, record_history)
% BLOCH_EVOLVE_ACTIVE_INDEXED 
% 利用预置的 .is_active 属性实现极速 Bloch 演化
% 
% --- 整合后的输入变量 ---
% Spins  - 结构体，包含：
%          .M      : [3 x N] 初始磁化矢量
%          .Pos    : [3 x N] 空间位置 (m)
%          .T1/T2  : [1 x N] 弛豫时间 (s)
%          .PD     : [1 x N] 质子密度 (M0)
%          .dB     : [1 x N] 静磁场不均匀性 (Hz)
% RF     - 结构体，包含：
%          .B1     : [1 x T] 复数射频序列 (T)
%          .df     : [1 x T] 射频发射偏置 (Hz)
%          .is_active :[1 x T]
% Grads  - 结构体，包含
%          .Gxyz [3 x T] 梯度场序列 (T/m)
%          .is_active :[1 x T]
% dt     - 步长 (s)

    % --- 1. 基本参数准备 ---
    gamma = 42.57e6;
    gamma_rad = 2 * pi * 42.57e6;
    [~, num_spins] = size(Spins.M);
    num_time = length(RF.B1);
    
    M = Spins.M;
    M_history = [];
    if record_history
        M_history = zeros(3, num_spins, num_time); 
    end

    % 预计算弛豫算子
    E1 = exp(-dt ./ Spins.T1);
    E2 = exp(-dt ./ Spins.T2);
    rec_PD = Spins.PD .* (1 - E1);
    
    % 预计算静磁场相位偏移 (phi_dB: 1 x N)
    phi_dB = 2 * pi * Spins.dB * dt;

    % --- 2. 核心演化循环 ---
    for t = 1:num_time
        
        % --- 获取当前时刻的主动标记 ---
        % RF.is_active(t) 和 Grads.is_active(t) 是预先定义好的布尔数组
        rf_on   = RF.is_active(t);
        grad_on = Grads.is_active(t);
        spatial_df = zeros(1,num_spins);
        if grad_on
                spatial_df = gamma * Grads.Gxyz(:,t)' * Spins.Pos;
        end
        % --- 分支路径选择 ---
        if rf_on
            % 路径 A: 射频活动期 (执行全 Rodrigues 旋转)
            % 这是唯一需要处理 3D 旋转和 RF 频率偏移 (df) 的地方
            df_total = Spins.dB - RF.df(t);
            if grad_on
                df_total = df_total + spatial_df; 
            end
            
            oz = 2 * pi * df_total;
            Omega = [repmat(real(RF.B1(t)) * gamma_rad, 1, num_spins); ...
                     repmat(imag(RF.B1(t)) * gamma_rad, 1, num_spins); ...
                     oz];
            M = apply_rodrigues(M, Omega, dt);

            % add 立即旋回到中心频率观察系（关键）
            phi = 2 * pi * RF.df(t) * dt;
            M = apply_z_rotation(M, phi);


        elseif grad_on
            % 路径 B: 仅梯度活动期 (仅 Z 轴平面旋转)
            % 绕过了复杂的 3D 旋转计算，直接更新相位
            phi_z = 2 * pi * (Spins.dB + spatial_df) * dt;
            M = apply_z_rotation(M, phi_z);

        else
            % 路径 C: 纯自由进动期 (Null 状态)
            % 仅受 B0 不均匀性影响，计算量极小
            M = apply_z_rotation(M, phi_dB);
        end

        % --- 3. 统一执行弛豫 ---
        M(1:2, :) = M(1:2, :) .* E2;
        M(3, :)   = M(3, :) .* E1 + rec_PD;

        if record_history
            M_history(:,:,t) = M; 
        end
    end
    M_final = M;
end

% --- 辅助函数：绕 Z 轴旋转 ---
function M = apply_z_rotation(M, phi)
    % phi 可以是标量或 [1 x N] 向量
    c = cos(phi); s = sin(phi);
    mx = M(1,:); my = M(2,:);
    M(1,:) = mx.*c - my.*s;
    M(2,:) = mx.*s + my.*c;
end

% --- 辅助函数：Rodrigues 3D 旋转 ---
function M = apply_rodrigues(M, Omega, dt)
    theta_v = sqrt(sum(Omega.^2, 1));
    theta = theta_v * dt;
    n = Omega ./ (theta_v + eps);
    c = cos(theta); s = sin(theta);
    dot_nM = sum(n .* M, 1);
    cross_nM = [n(2,:).*M(3,:) - n(3,:).*M(2,:); ...
                n(3,:).*M(1,:) - n(1,:).*M(3,:); ...
                n(1,:).*M(2,:) - n(2,:).*M(1,:)];
    M = M.*c + cross_nM.*s + n.*(dot_nM.*(1-c));
end