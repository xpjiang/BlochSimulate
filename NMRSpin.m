classdef NMRSpin < handle
    properties
        PD          % 质子密度
        T1          % s
        T2          % s
        dB          % 局部磁场偏差 (Hz)
        Position    % [x, y, z] (m)
        M           % [Mx; My; Mz]
    end
    
    methods
        function obj = NMRSpin(pd, t1, t2, db, pos)
            obj.PD = pd;
            obj.T1 = t1;
            obj.T2 = t2;
            obj.dB = db;
            obj.Position = pos;
            obj.M = [0; 0; pd]; 
        end
        
        % --- 核心演化函数 (The Evolve Function) ---
        % --- 考虑偏共振的演化函数 ---
        function evolve(obj, dt, B1, Grads, df_tx)
            % dt: 时间步长 (s)
            % B1: 复数射频场 (Tesla)
            % Grads: 当前梯度场 [Gx, Gy, Gz] (T/m)
            % df_tx: 射频发射中心频率偏移 (Hz) - 新增参数
            
            gamma = 42.57e6;
            gamma_rad = gamma * 2 * pi;
            
            % 1. 射频场分量 (B1 作用在横平面)
            omega_x = gamma_rad * real(B1);
            omega_y = gamma_rad * imag(B1);
            
            % 2. 纵向有效频率偏差 (核心修改)
            % 自旋物理频率 f_spin = dB + f_grad
            % 在射频参考系下的频率 = f_spin - df_tx
            f_grad = gamma * sum(obj.Position .* Grads);
            f_eff = (obj.dB + f_grad) - df_tx; 
            omega_z = 2 * pi * f_eff;
            
            % 3. 构建布洛赫矩阵 (4x4 增广矩阵法)
            R1 = 1/obj.T1;
            R2 = 1/obj.T2;
            
            % 注意：omega_z 现在包含了所有的偏共振效应
            A = [ -R2,      omega_z,  -omega_y;
                  -omega_z,  -R2,       omega_x;
                   omega_y,  -omega_x, -R1 ];
            
            C = [0; 0; obj.PD/obj.T1];
            F = [A, C; 0 0 0 0];
            
            % 4. 演化更新
            ExpFdt = expm(F * dt);
            M_aug = ExpFdt * [obj.M; 1];
            obj.M = M_aug(1:3);
        end
    end
end