% question2_x1.m
% Generate x1[n], compute autocorrelation, correlogram and compare periodogram
close all; clear; clc;

N = 1024;
sigma1 = sqrt(1/26);
w1 = sigma1 * randn(1, N);

% Define the coefficients and create x1
coefficients = [1 -3 -4];
x1 = filter(coefficients, 1, w1);

% Autocorrelation (biased)
r = xcorr(x1, 'biased');
lags = -N+1 : N-1;

% Plot autocorrelation
figure;
plot(lags, r);
title('Biased autocorrelation of x1[n]');
xlabel('Lag');
ylabel('Autocorrelation');
grid on;

% Correlogram (compute PSD from biased autocorrelation)
omega = linspace(0, pi, 2049);
Pxx = zeros(1, length(omega));
for k = 1:length(omega)
    Pxx(k) = sum(r .* exp(-1j * omega(k) * lags));
end

% Plot the correlogram
figure;
plot(omega, abs(Pxx));
title('Correlogram of x1[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power Spectral Density');
xlim([0 pi]);
grid on;

% Compare theoretical spectrum, periodogram and correlogram on one figure
w = linspace(0, pi, 2049);
sxx1 = (26 - 8*cos(2*w) + 18*cos(w)) / 26;

figure;
hold on;
plot(w, sxx1, 'k', 'DisplayName', 'Theoretical Spectrum of X_1[n]');

% Periodogram
N2 = 4096;
x1p = [x1, zeros(1, N2 - N)];
X1 = fft(x1p);
Pxx_per = (1/N) * abs(X1).^2;
freq = linspace(0, 2*pi, N2);
plot(freq, Pxx_per, 'b', 'DisplayName', 'Periodogram of X_1[n]');

% Correlogram (recompute with fresh realization)
w1 = sigma1 * randn(1, N);
x1 = filter(coefficients, 1, w1);
r = xcorr(x1, 'biased');
Pxx_corr = zeros(1, length(omega));
for k = 1:length(omega)
    Pxx_corr(k) = sum(r .* exp(-1j * omega(k) * lags));
end
plot(omega, abs(Pxx_corr), 'g', 'DisplayName', 'Correlogram of X_1[n]');

title('Spectrum, Periodogram and Correlogram of X_1[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
legend;
grid on;
hold off;
