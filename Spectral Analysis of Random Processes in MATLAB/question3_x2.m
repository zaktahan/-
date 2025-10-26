% question3_x2.m
% Monte-Carlo average of spectral estimators for X2 (AR process)
close all; clear; clc;

% Theoretical spectrum
w = linspace(0, pi, 2049);
sxx2 = 0.51 ./ (1.49 - 1.4 * cos(w));
sxx2 = abs(sxx2);

% Parameters
N = 1024;
num_Monte_Carlo = 100;
M = 4096;
f_samples = M/2 + 1;
sigma2 = sqrt(0.51);

% Preallocations
Pxx_A_all = zeros(f_samples, num_Monte_Carlo);
Pxx_B_all = zeros(f_samples, num_Monte_Carlo);
Pxx_welch_all = zeros(f_samples, num_Monte_Carlo);
Pxx_welch_all_B = zeros(f_samples, num_Monte_Carlo);

% Bartlett L=64, K=16 and Welch L=64 variant
L_A = 64; K_A = 16;
L_welch = 64; D_welch = 48; K_welch = 61;

for r = 1:num_Monte_Carlo
    w2 = sigma2 * randn(N,1);
    X2 = zeros(N,1);
    for n = 2:N
        X2(n) = 0.7 * X2(n-1) + w2(n);
    end
    Pxx_A_all(:, r) = bartlett_method(X2, L_A, K_A, M);

    x2 = zeros(N,1);
    for n = 2:N
        x2(n) = 0.7 * x2(n-1) + w2(n);
    end
    Pxx_welch_all(:, r) = welch_method(x2, L_welch, D_welch, K_welch, M);
end

% Bartlett L=16, K=64 and Welch L=16 variant
L_B = 16; K_B = 64;
L_welch_B = 16; D_welch_B = 12; K_welch_B = 253;

for r = 1:num_Monte_Carlo
    w2 = sigma2 * randn(N,1);
    X2 = zeros(N,1);
    for n = 2:N
        X2(n) = 0.7 * X2(n-1) + w2(n);
    end
    Pxx_B_all(:, r) = bartlett_method(X2, L_B, K_B, M);

    x2 = zeros(N,1);
    for n = 2:N
        x2(n) = 0.7 * x2(n-1) + w2(n);
    end
    Pxx_welch_all_B(:, r) = welch_method(x2, L_welch_B, D_welch_B, K_welch_B, M);
end

meanPxx_A = mean(Pxx_A_all, 2);
meanPxx_B = mean(Pxx_B_all, 2);
meanPxx_welch = mean(Pxx_welch_all, 2);
meanPxx_welch_B = mean(Pxx_welch_all_B, 2);

% Blackman-Tukey small-l estimates (l=2 and l=4)
N2 = 4098;
Mc = 100;
sigma2 = sqrt(0.51);
Pxx_mean = zeros(1, N2);
Pxx_mean_4 = zeros(1, N2);

for m = 1:Mc
    w2 = sigma2 * randn(1, N);
    x2 = zeros(1, N);
    a = 0.7;
    for n = 2:N
        x2(n) = a * x2(n-1) + w2(n);
    end
    r = xcorr(x2, 'biased');
    r = circshift(r, -1023 + 2);
    r = r(1:2*2+1);
    numElementsToPad = N2 - length(r);
    r = [r, zeros(1, numElementsToPad)];
    s = fft(r);
    Pxx_mean = Pxx_mean + s;
end
Pxx_mean = Pxx_mean / Mc;

for m = 1:Mc
    w2 = sigma2 * randn(1, N);
    x2 = zeros(1, N);
    a = 0.7;
    for n = 2:N
        x2(n) = a * x2(n-1) + w2(n);
    end
    r_4 = xcorr(x2, 'biased');
    r_4 = circshift(r_4, -1023 + 4);
    r_4 = r_4(1:2*4+1);
    numElementsToPad = N2 - length(r_4);
    r_4 = [r_4, zeros(1, numElementsToPad)];
    s_4 = fft(r_4);
    Pxx_mean_4 = Pxx_mean_4 + s_4;
end
Pxx_mean_4 = Pxx_mean_4 / Mc;

% Frequency vectors and plot
f_samples = M/2 + 1;
omega = linspace(0, 2*pi, N2);
f = (0:f_samples-1) * (pi / (f_samples-1));

figure;
plot(f, meanPxx_A, 'DisplayName', 'Bartlett L=64');
hold on;
plot(w, sxx2, 'DisplayName', 'Theoretical Spectrum X_2[n]');
plot(f, meanPxx_B, 'DisplayName', 'Bartlett L=16');
plot(f, meanPxx_welch, 'DisplayName', 'Welch L=64, D=48');
plot(f, meanPxx_welch_B, 'DisplayName', 'Welch L=16, D=12');
plot(omega, abs(Pxx_mean), 'DisplayName', 'BT of X_2[n] l = 2');
plot(omega, abs(Pxx_mean_4), 'DisplayName', 'BT of X_2[n] l = 4');
hold off;
title('Average of the spectrum value for X_2[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
legend('show');
grid on;

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
    N = length(x);
    Pxx = zeros(M/2 + 1, 1);
    for k = 1:K
        start_idx = (k-1) * (L - D) + 1;
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
