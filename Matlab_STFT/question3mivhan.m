
clc; clear; close all;

% Section A: Define Parameters and Analysis Window
M = 16; % Number of frequency bands
N = 10*M; % Length of analysis window
n = -N:N; % Time index
h = sin(pi * n / M) ./ (pi * n); 
h(N+1) = 1/M;

%% Section A
W = linspace(-pi, pi, 1024); % Frequency range
H_dtft = sum(h.' .* exp(-1j * (n.') * W), 1); % Compute DTFT
figure;
plot(W, abs(H_dtft));
xlabel('\omega'); ylabel('|H(e^j\omega)|');
title('DTFT of Analysis Window h[n]');
grid on;
figure;
hold on;
for k0 = 0:M-1
    omega_shift = W - (2*pi*k0/M); % Shift frequency by 2*pi*k0/M
    H_k0 = interp1(W, H_dtft, omega_shift, 'linear', 0); % Interpolate shifted DTFT
    plot(W, abs(H_k0), 'DisplayName', sprintf('k0 = %d', k0));
end
xlabel('\omega'); ylabel('|H_{k0}(e^j\omega)|');
title('QMF H_{k0} Response for All k0');
legend;
grid on;
hold off;

%% Section D: Generate Signal x[n] = cos(2*pi*k*n/M)
k = 4; % Example k value from exam
x = cos(2*pi*k*n/M);

% Section D: Compute and Plot Spectrogram
window = hamming(256);
overlap = 128;
nfft = 512;
figure;
spectrogram(x, window, overlap, nfft, M, 'yaxis');
title('Spectrogram of x[n]');
colorbar;
