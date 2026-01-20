function ax = plot_mri_sequence2(Spins, RF, Grads, t_ms, signal)
%PLOT_MRI_SEQUENCE  改进版：支持多轴梯度显示与 ADC 采样窗标注
%
% 输入：
%   Spins.Pos   [3 x N] (m)
%   Spins.dB    [1 x N] (Hz)
%   RF.B1       [1 x T] (T, complex)
%   RF.is_adc   [1 x T] (Logical, ADC 采样窗标记) <- 新增建议输入
%   Grads.Gxyz  [3 x T] (T/m)
%   t_ms        [1 x T] 时间 (ms)
%   signal      [1 x T] 复信号

    figure('Color', 'w');
    % --- 紧凑布局 ---
    tiledlayout(4,1,'TileSpacing','compact','Padding','compact')

    %% --- 子图 1：自旋分布 (保持不变) ---
    ax.ax1 = nexttile;
    scatter3(Spins.Pos(1,:)*1e3, Spins.Pos(2,:)*1e3, Spins.Pos(3,:)*1e3, ...
             15, abs(Spins.dB), 'filled');
    axis equal; grid on;
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
    title('Spins Distribution (|dB| colored)');
    colormap(gca, 'jet'); colorbar;

    %% --- 子图 2：RF + ADC Window ---
    ax.ax2 = nexttile;
    hold on;
    % 绘制 ADC 采样窗背景 (如果有提供 RF.is_adc)
    if isfield(RF, 'is_adc') && any(RF.is_adc)
        yl = [-1, 1] * max(abs(RF.B1)*1e6) * 1.2; % 动态获取高度
        if all(yl == 0), yl = [0 1]; end
        % 寻找连贯的 ADC 区间
        adc_starts = find(diff([0, RF.is_adc]) == 1);
        adc_ends = find(diff([RF.is_adc, 0]) == -1);
        for i = 1:length(adc_starts)
            patch([t_ms(adc_starts(i)), t_ms(adc_ends(i)), t_ms(adc_ends(i)), t_ms(adc_starts(i))], ...
                  [yl(1) yl(1) yl(2) yl(2)], [1 1 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'on', 'DisplayName', 'ADC Window');
        end
    end
    plot(t_ms, abs(RF.B1)*1e6, 'k', 'LineWidth', 1.2, 'DisplayName', '|B1|');
    plot(t_ms, real(RF.B1)*1e6, 'LineWidth', 1.0, 'DisplayName', 'Real(B1)');
    plot(t_ms, imag(RF.B1)*1e6, 'LineWidth', 1.0, 'DisplayName', 'Imag(B1)');
    ylabel('RF (\muT)');
    title('RF Pulses & ADC Window');
    grid on; legend('Location', 'northeast');
    hold off;

    %% --- 子图 3：G_x, G_y, G_z 多轴梯度 ---
    ax.ax3 = nexttile;
    plot(t_ms, Grads.Gxyz(1,:)*1e3, 'r', 'LineWidth', 1.2, 'DisplayName', 'G_x'); hold on;
    plot(t_ms, Grads.Gxyz(2,:)*1e3, 'g', 'LineWidth', 1.2, 'DisplayName', 'G_y');
    plot(t_ms, Grads.Gxyz(3,:)*1e3, 'b', 'LineWidth', 1.2, 'DisplayName', 'G_z');
    ylabel('G (mT/m)');
    title('Gradients (X, Y, Z)');
    grid on; legend('Location', 'northeast');
    hold off;

    %% --- 子图 4：Signal ---
    ax.ax4 = nexttile;
    plot(t_ms, abs(signal), 'k', 'LineWidth', 1.5, 'DisplayName', 'Magnitude'); hold on;
    plot(t_ms, real(signal), 'b', 'LineWidth', 1.0, 'DisplayName', 'Real');
    plot(t_ms, imag(signal), 'r', 'LineWidth', 1.0, 'DisplayName', 'Imag');
    ylabel('Signal (M_{xy})');
    xlabel('Time (ms)');
    title('Transverse Magnetization Evolution');
    grid on; legend('Location', 'northeast');
    hold off;

    %% --- 共享 x 轴与格式化 ---
    linkaxes([ax.ax2, ax.ax3, ax.ax4], 'x');
    set([ax.ax2, ax.ax3], 'XTickLabel', [], 'XColor', 'none');
    
    % 统一留白
    expand_axes_y([ax.ax2, ax.ax3, ax.ax4], 0.1);
end

function expand_axes_y(axs, ratio)
    for k = 1:numel(axs)
        yl = ylim(axs(k));
        dy = diff(yl);
        if dy == 0, dy = 1; end
        ylim(axs(k), yl + ratio*[-dy dy]);
    end
end