clc; clear; close all;

%% Section A: Define Parameters and Analysis Window
M = 16;                  % Number of frequency bands
N = 10*M;                % Half-length of analysis window
n = -N:N;                % Time indices
h = sin(pi * n / M) ./ (pi * n); 
h(N+1) = 1/M;            % Handle division by zero at n=0

%% Section B: Compute DTFT of Analysis Window
W = linspace(-pi, pi, 1024);          % Frequency axis
H_dtft = sum(h.' .* exp(-1j * (n.') * W), 1);  % Compute DTFT

% Plot magnitude of DTFT
figure;
plot(W, abs(H_dtft), 'LineWidth', 1.5);
xlabel('\omega'); ylabel('|H(e^{j\omega})|');
title('DTFT of Analysis Window h[n]');
grid on;

%% Section C: QMF Responses for All k0
figure; hold on;
for k0 = 0:M-1
    omega_shift = W - (2*pi*k0/M);                   % Shift frequency by 2*pi*k0/M
    H_k0 = interp1(W, H_dtft, omega_shift, 'linear', 0);  % Interpolate shifted DTFT
    plot(W, abs(H_k0), 'DisplayName', sprintf('k0 = %d', k0));
end
xlabel('\omega'); ylabel('|H_{k0}(e^{j\omega})|');
title('QMF H_{k0} Response for All k0');
legend show; grid on; hold off;

%% Section D: Generate Signal x[n] = cos(2*pi*k*n/M)
k = 4;                  % Example k value
x = cos(2*pi*k*n/M);

%% Section E: Compute and Plot Spectrogram
window = hamming(256);  % Analysis window for spectrogram
overlap = 128;          % Number of overlapping samples
nfft = 512;             % FFT points
figure;
spectrogram(x, window, overlap, nfft, M, 'yaxis');
title('Spectrogram of x[n]');
colorbar;
