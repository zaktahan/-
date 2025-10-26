% question3_x1.m
% Monte-Carlo average of spectral estimators for X1 (MA process)
close all; clear; clc;

% Theoretical spectrum
w = linspace(0, pi, 2049);
sxx1 = (26 - 8*cos(2*w) + 18*cos(w)) / 26;

% Parameters
N = 1024;
num_Monte_Carlo = 100;
M = 4096;
f_samples = M/2 + 1;
sigma1 = sqrt(1/26);

% Preallocate results
Pxx_A_all = zeros(f_samples, num_Monte_Carlo);
Pxx_B_all = zeros(f_samples, num_Monte_Carlo);
Pxx_welch_all = zeros(f_samples, num_Monte_Carlo);
Pxx_welch_additional_all = zeros(f_samples, num_Monte_Carlo);

% Bartlett L=64, K=16
L_A = 64; K_A = 16;
for r = 1:num_Monte_Carlo
    w1 = sigma1 * randn(N,1);
    x1 = zeros(N,1);
    for n = 3:N
        x1(n) = w1(n) - 3*w1(n-1) - 4*w1(n-2);
    end
    Pxx_A_all(:, r) = bartlett_method(x1, L_A, K_A, M);
end
meanPxx_A = mean(Pxx_A_all, 2);

% Bartlett L=16, K=64
L_B = 16; K_B = 64;
for r = 1:num_Monte_Carlo
    w1 = sigma1 * randn(N,1);
    x1 = zeros(N,1);
    for n = 3:N
        x1(n) = w1(n) - 3*w1(n-1) - 4*w1(n-2);
    end
    Pxx_B_all(:, r) = bartlett_method(x1, L_B, K_B, M);
end
meanPxx_B = mean(Pxx_B_all, 2);

% Welch variations
L_welch = 64; D_welch = 48; K_welch = 61;
for r = 1:num_Monte_Carlo
    w1 = sigma1 * randn(N,1);
    x1 = zeros(N,1);
    for n = 3:N
        x1(n) = w1(n) - 3*w1(n-1) - 4*w1(n-2);
    end
    Pxx_welch_all(:, r) = welch_method(x1, L_welch, D_welch, K_welch, M);
end
meanPxx_welch = mean(Pxx_welch_all, 2);

L_welch_add = 16; D_welch_add = 12; K_welch_add = 253;
for r = 1:num_Monte_Carlo
    w1 = sigma1 * randn(N,1);
    x1 = zeros(N,1);
    for n = 3:N
        x1(n) = w1(n) - 3*w1(n-1) - 4*w1(n-2);
    end
    Pxx_welch_additional_all(:, r) = welch_method(x1, L_welch_add, D_welch_add, K_welch_add, M);
end
meanPxx_welch_additional = mean(Pxx_welch_additional_all, 2);

% Blackman-Tukey style (small-l)
N2 = 4098;
l_values = [2, 4];
Mc = 100;
Pxx_mean = zeros(2, 2049);
for idx = 1:2
    l = l_values(idx);
    for m = 1:Mc
        w1 = sigma1 * randn(1, N);
        coefficients = [1 -3 -4];
        x1 = filter(coefficients, 1, w1);
        r = xcorr(x1, 'biased');
        r = circshift(r, -1023 + l);
        r = r(1:2*l+1);
        numElementsToPad = N2 - length(r);
        r = [r zeros(1, numElementsToPad)];
        s = fft(r);
        s = s(1:2049);
        Pxx_mean(idx, :) = Pxx_mean(idx, :) + s;
    end
    Pxx_mean(idx, :) = Pxx_mean(idx, :) / Mc;
end

omega = linspace(0, pi, 2049);

% Plot averages
figure; hold on;
plot(w, sxx1, 'DisplayName', 'Theoretical Spectrum X_1[n]');
plot((0:f_samples-1) * (pi/(f_samples-1)), meanPxx_A, 'DisplayName', 'Bartlett L = 64');
plot((0:f_samples-1) * (pi/(f_samples-1)), meanPxx_B, 'DisplayName', 'Bartlett L = 16');
plot((0:f_samples-1) * (pi/(f_samples-1)), meanPxx_welch, 'DisplayName', 'Welch L = 64, D=48');
plot((0:f_samples-1) * (pi/(f_samples-1)), meanPxx_welch_additional, 'DisplayName', 'Welch L=16, D=12');
plot(omega, abs(Pxx_mean(2, :)), 'DisplayName', 'BT of X_1[n] l = 4');
plot(omega, abs(Pxx_mean(1, :)), 'DisplayName', 'BT of X_1[n] l = 2');
title('Average of the spectrum value for X_1[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
legend;
grid on;
hold off;

% --- Auxiliary functions ---
function [Pxx] = bartlett_method(x, L, K, M)
    Pxx = zeros(M/2 + 1, 1);
    for k = 1:K
        segment = x((k-1)*L + 1 : k*L);
        segment_fft = fft(segment, M);
        periodogram = (1 / L) * abs(segment_fft(1:M/2 + 1)).^2;
        Pxx = Pxx + periodogram;
    end
    Pxx = Pxx / K;
end

function [Pxx] = welch_method(x, L, D, K, M)
    Pxx = zeros(M/2 + 1, 1);
    N = length(x);
    for k = 1:K
        start_idx = (k-1)*(L-D) + 1;
        end_idx = start_idx + L - 1;
        if end_idx > N
            break;
        end
        segment = x(start_idx:end_idx);
        windowed_segment = segment;
        segment_fft = fft(windowed_segment, M);
        periodogram = (1 / L) * abs(segment_fft(1:M/2 + 1)).^2;
        Pxx = Pxx + periodogram;
    end
    Pxx = Pxx / k;
end
