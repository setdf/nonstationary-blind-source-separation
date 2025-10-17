clear; clc; close all;
rng(3)
c = [0.2 0.4 0.6 -0.1 -0.3];
d = [0.1 0.3 -0.2 0.5 -0.3];
fs = 20;
A = [0.8 -0.6; 0.6 0.8];
num_of_windows = 5;


s1 = [];
s2 = [];
for k = 1 : num_of_windows
    t = (k - 1) : 1 / fs : k - (1 / fs);
    s1 = [s1 c(k) * sin(2 * pi * t)];
    s2 = [s2 d(k) * sin(4 * pi * t)];
end
t_total = 0 : 1 / fs : num_of_windows - (1 / fs);
S = [s1; s2];
X = A * S;

figure
plot(t_total, S(1, :), 'DisplayName', 's_1')
hold on
plot(t_total, S(2, :), 'DisplayName', 's_2')
hold off
legend
xlabel('time(s)')
title('Sources')

figure
plot(t_total, X(1, :), 'DisplayName', 'x_1')
hold on
plot(t_total, X(2, :), 'DisplayName', 'x_2')
hold off
legend
xlabel('time(s)')
title('X')

%% b
window_size = size(X, 2) / num_of_windows;
Rxs = compute_Rxs(X, num_of_windows, window_size);
two_windows(Rxs(1:2), X, S)

%% c
[U, D] = eig(X * X');
Z = D ^ -0.5 * U' * X;
Rzs = compute_Rxs(Z, num_of_windows, window_size);
error = all_windows(Rzs, S, Z);
disp(['Error using all windows: ', num2str(error)]);
%% d
disp('Noisy data:')
SNR_desired_linear = 100;
Y = generate_noisy_signal(X, SNR_desired_linear);
Rys = compute_Rxs(Y, num_of_windows, window_size);
two_windows(Rys, Y, S)
[U, D] = eig(Y * Y');
Z = D ^ -0.5 * U' * Y;
Rzs = compute_Rxs(Z, num_of_windows, window_size);
error = all_windows(Rzs, S, Z);
disp(['Error using all windows: ', num2str(error)]);
%% e
errors = zeros(1, 4);
for k = 2:5
    error_sum = 0;
    for trial = 1:100
        Y = generate_noisy_signal(X, SNR_desired_linear);
        [U, D] = eig(Y * Y');
        Z = D ^ -0.5 * U' * Y;
        Rzs = compute_Rxs(Z, num_of_windows, window_size);
        error = all_windows(Rzs(1:k), S, Z); 
        error_sum = error_sum + error;
    end
    errors(k-1) = error_sum / 100;
end

figure;
plot(2:5, errors, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('i (number of windows)');
ylabel('mean error');
title('mean errors for different number of windows');
%% f
SNR_desired_db = [5 10 15 20];
SNRs_desired_linear = 10.^(SNR_desired_db / 10);

errors = zeros(1, 4);
for i = 1:4
    error_sum = 0;
    for trial = 1:100
        Y = generate_noisy_signal(X, SNRs_desired_linear(i));
        [U, D] = eig(Y * Y');
        Z = D ^ -0.5 * U' * Y;
        Rzs = compute_Rxs(Z, num_of_windows, window_size);
        error = all_windows(Rzs, S, Z); 
        error_sum = error_sum + error;
    end
    errors(i) = error_sum / 100;
end

figure;
plot(SNR_desired_db, errors, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)');
ylabel('mean error');
title('mean errors for different SNR levels');


%% Functions
function Y = generate_noisy_signal(X, SNR_desired_linear)
    noise = randn(size(X));
    W = noise / norm(noise, 'fro');

    signal_power = norm(X, 'fro')^2;
    noise_power = signal_power / SNR_desired_linear;
    sigma = sqrt(noise_power);

    Y = X + sigma * W;
end

function Rxs = compute_Rxs(X, num_of_windows, window_size)
    Rxs = cell(1, num_of_windows);
    for k = 1:num_of_windows
        idx_start = round((k - 1) * window_size) + 1;
        idx_end   = round(k * window_size);
        Xk = X(:, idx_start:idx_end);
        Rxs{k} = (Xk * Xk');
    end
end

function result = R(B, Rx_all, m)
    [n, ~] = size(B);
    result = 0;
    for k = 1:length(Rx_all)
        Rk = Rx_all{k};
        for j = 1:n
            if j == m
                continue
            end
            bj = B(:, j);
            temp = Rk * bj;
            result = result + temp * temp';
        end
    end
end

function B = optimize_B(Rxs)
    B_initial = ones(2, 2);
    cost_old = [Inf, Inf];
    B = B_initial;
    optimized = false;    
    while ~optimized
        for m = 1 : 2
            R_m = R(B, Rxs, m);
            [b, D] = eig(R_m);
            [~, min_idx] = min(diag(D));
            b_m = b(:, min_idx);
            if m == 2
                b1 = B(:, 1);
                b_m = (eye(size(b1, 1)) - b1 * b1') * b_m;
                b_m = b_m / norm(b_m);
            end
            cost_new(m) = b_m' * R_m * b_m;
            if cost_new(m) >= cost_old(m)
                optimized = true;
            end
            B(:, m) = b_m;
            cost_old(m) = cost_new(m);
        end
    end
end

function two_windows(Rxs, X, S)
    [Bt, L] = eig(Rxs{1}, Rxs{2});
    S_estimated = Bt' * X;

    S_norm = (S - mean(S, 2)) ./ std(S, 0, 2);
    S_est_norm = -(S_estimated - mean(S_estimated, 2)) ./ std(S_estimated, 0, 2);
    S_est_flipped = S_est_norm([2 1], :);

    error = norm(S_est_flipped - S_norm, 'fro')^2 / norm(S_norm, 'fro')^2;
    disp(['Error using two windows: ', num2str(error)]);
end

function error = all_windows(Rxs, S, X)
    B = optimize_B(Rxs);

    S_estimated = B' * X;

    S_norm = (S - mean(S, 2)) ./ std(S, 0, 2);
    S_est_norm = (S_estimated - mean(S_estimated, 2)) ./ std(S_estimated, 0, 2);
    S_est_flipped = S_est_norm([2 1], :);

    error_original = norm(S_est_norm - S_norm, 'fro')^2 / norm(S_norm, 'fro')^2;
    error_flipped = norm(S_est_norm([2 1], :) - S_norm, 'fro')^2 / norm(S_norm, 'fro')^2;
    
    if error_flipped < error_original
        error = error_flipped;
    else
        error = error_original;
    end
    
end


