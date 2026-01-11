function ax = plot_mri_sequence(Spins, RF, Grads, t_ms, signal)
%PLOT_MRI_SEQUENCE  MRI 序列与信号可视化（一行调用）
%
% ax = plot_mri_sequence(Spins, RF, Grads, t_ms, signal)
%
% 输入：
%   Spins.Pos   [3 x N]  (m)
%   Spins.dB    [1 x N]  (Hz)
%   RF.B1       [1 x T]  (T, complex)
%   Grads.Gxyz  [3 x T]  (T/m)
%   t_ms        [1 x T]  时间 (ms)
%   signal      [1 x T]  复信号
%
% 输出：
%   ax          axes 句柄结构体

    figure

    % --- 紧凑布局（关键） ---
    tiledlayout(4,1,'TileSpacing','compact','Padding','compact')

    %% --- 子图 1：自旋分布 ---
    ax.ax1 = nexttile;
    scatter3(Spins.Pos(1,:)*1e3, ...
             Spins.Pos(2,:)*1e3, ...
             Spins.Pos(3,:)*1e3, ...
             15, abs(Spins.dB), 'filled');
    axis equal
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
    title('Spins Distribution (|dB| colored)')
    colormap jet
    colorbar

    %% --- 子图 2：RF ---
    ax.ax2 = nexttile;
    plot(t_ms, abs(RF.B1)*1e6, 'LineWidth', 1.2); hold on
    plot(t_ms, real(RF.B1)*1e6, 'LineWidth', 1.2)
    plot(t_ms, imag(RF.B1)*1e6, 'LineWidth', 1.2)
    ylabel('RF (\muT)')
    grid on
    hold off;

    %% --- 子图 3：Gradient ---
    ax.ax3 = nexttile;
    plot(t_ms, Grads.Gxyz(1,:)*1e3, 'LineWidth', 1.2)
    ylabel('G_x (mT/m)')
    grid on

    %% --- 子图 4：Signal ---
    ax.ax4 = nexttile;
    plot(t_ms, abs(signal), 'LineWidth', 1.5); hold on
    plot(t_ms, real(signal), 'LineWidth', 1.2)
    plot(t_ms, imag(signal), 'LineWidth', 1.2)
    ylabel('M_{xy}')
    xlabel('Time (ms)')
    grid on;
	hold on;

    %% --- 共享 x 轴（2–4） ---
    linkaxes([ax.ax2, ax.ax3, ax.ax4], 'x')

    % --- 只显示最下面的 x 轴 ---
    set([ax.ax2, ax.ax3], 'XTickLabel', [], 'XColor', 'none')

    %% --- 统一留白（y 轴） ---
    expand_axes_y([ax.ax2, ax.ax3, ax.ax4], 0.03)
end


%% ===== 辅助函数：y 轴自动留白 =====
function expand_axes_y(axs, ratio)
    for k = 1:numel(axs)
        yl = ylim(axs(k));
        dy = diff(yl);
        ylim(axs(k), yl + ratio*[-dy dy]);
    end
end
