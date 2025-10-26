% question2_x2.m
% Generate AR(1) x2[n], compute autocorrelation, correlogram and compare periodogram
close all; clear; clc;

N = 1024;
sigma2 = sqrt(0.51);
w2 = sigma2 * randn(1, N);

% Generate AR(1) process x2(n) = 0.7 x2(n-1) + w2(n)
x2 = zeros(1, N);
a = 0.7;
for n = 2:N
    x2(n) = a * x2(n-1) + w2(n);
end

% Autocorrelation (biased)
r = xcorr(x2, 'biased');
lags = -N+1 : N-1;

% Plot autocorrelation
figure;
plot(lags, abs(r));
title('Biased autocorrelation of x2[n]');
xlabel('Lag');
ylabel('Autocorrelation (abs)');
grid on;

% Correlogram
omega = linspace(0, pi, 2049);
Pxx = zeros(1, length(omega));
for k = 1:length(omega)
    Pxx(k) = sum(r .* exp(-1j * omega(k) * lags));
end

figure;
plot(omega, abs(Pxx));
title('Correlogram of x2[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power Spectral Density');
xlim([0 pi]);
grid on;

% Theoretical spectrum
w = linspace(0, pi, 2049);
sxx2 = 0.51 ./ (1.49 - 1.4 * cos(w));

% Periodogram comparison
figure;
hold on;
plot(w, sxx2, 'g', 'DisplayName', 'Theoretical Spectrum of X_2[n]');

N2 = 4096;
x2p = [x2, zeros(1, N2 - N)];
X2 = fft(x2p);
Pxx_per = (1/N) * abs(X2).^2;
freq = linspace(0, 2*pi, N2);
plot(freq, Pxx_per, 'b', 'DisplayName', 'Periodogram of X_2[n]');

% Correlogram (another realization)
w2 = sigma2 * randn(1, N);
x2 = zeros(1, N);
for n = 2:N
    x2(n) = a * x2(n-1) + w2(n);
end
r = xcorr(x2, 'biased');
Pxx_corr = zeros(1, length(omega));
for k = 1:length(omega)
    Pxx_corr(k) = sum(r .* exp(-1j * omega(k) * lags));
end
plot(omega, abs(Pxx_corr), 'b', 'DisplayName', 'Correlogram of X_2[n]');

title('Spectrum, Periodogram and Correlogram of X_2[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
legend;
grid on;
hold off;
