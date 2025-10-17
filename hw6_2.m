clear; clc; close all;
load hw6-X1.mat
load hw6-X2.mat

fs = 100;
num_of_windows = 2;

source_estimation(X1, num_of_windows, fs)
source_estimation(X2, num_of_windows, fs)

function source_estimation(X, num_of_windows, fs)
    tau = 5;
    Rxs = compute_Rxs(X, num_of_windows, 200, tau);

    [Bt, L] = eig(Rxs{1}, Rxs{2});
    S_estimated = Bt' * X;
    
    t_total = 0 : 1/fs : size(X, 2) / fs - 1/fs;

    figure;
    plot(t_total, S_estimated(1, :), 'DisplayName', '$\hat{s}_1$');
    hold on;
    plot(t_total, S_estimated(2, :), 'DisplayName', '$\hat{s}_2$');
    hold off;
    legend('Interpreter', 'latex');
    xlabel('Time (s)');
    title('Estimated sources');
    
    S_estimated_fft_1 = fft(S_estimated(1, :));
    S_estimated_fft_2 = fft(S_estimated(2, :));
    f = (0:length(S_estimated)-1) * fs / length(S_estimated);
    
    figure;
    subplot(2, 1, 1);
    plot(f(1:floor(length(f)/2)), abs(S_estimated_fft_1(1:floor(length(S_estimated_fft_1)/2))));
    title('fourier transform of $\hat{s}_1$', 'Interpreter', 'latex');
    xlabel('frequency(Hz)');
    
    subplot(2, 1, 2);
    plot(f(1:floor(length(f)/2)), abs(S_estimated_fft_2(1:floor(length(S_estimated_fft_2)/2))));
    title('fourier transform of $\hat{s}_2$', 'Interpreter', 'latex');
    xlabel('frequency(Hz)');
end

function Rxs = compute_Rxs(X, num_of_windows, window_size, tau)
    Rxs = cell(1, num_of_windows);
    first_window = X(:, 1:window_size);
    for k = 1:num_of_windows
        idx_start = round((k) * tau) + 1;
        idx_end = idx_start + window_size - 1;
        
        Xk = X(:, idx_start:idx_end);
        Rxs{k} = (first_window * Xk');
    end
end
