% question3_4_x1.m
% Mean periodogram and bias/variance/MSE calculations for X1
close all; clear; clc;

Mc = 100;
N = 1024;
N2 = 4096;
M = 4096;
sigma1 = sqrt(1/26);
Pxx_mean = zeros(1, N2);

for m = 1:Mc
    w1 = sigma1 * randn(1, N);
    coefficients = [1 -3 -4];
    x1 = filter(coefficients, 1, w1);
    x1 = [x1, zeros(1, N2 - N)];
    X1 = fft(x1);
    Pxx = (1/N) * abs(X1).^2;
    Pxx_mean = Pxx_mean + Pxx;
end
Pxx_mean = Pxx_mean / Mc;

freq = linspace(0, 2*pi, N2);

% Plot the mean periodogram
figure;
subplot(2,2,1);
plot(freq(1:N2/2), (Pxx_mean(1:N2/2)), 'LineWidth', 2);
title('Mean Periodogram Spectrum (Monte Carlo) - X1');
xlabel('Frequency (w)');
ylabel('Spectrum');
grid on;
xlim([0, pi]);

% Theoretical S
S = (26 - 8*cos(freq) + 9*cos(freq)) / 26; % kept original expression from code
B = Pxx_mean - S;

subplot(2,2,2);
plot(freq(1:N2/2), abs(B(1:N2/2)), 'LineWidth', 2);
title('Estimated - Theoretical (|B(e^{j\omega})|) - X1');
xlabel('Frequency (w)');
ylabel('|B(e^{j\omega})|');
grid on;
xlim([0, pi]);

% Variance
var = zeros(1, N2);
for m = 1:Mc
    w1 = sigma1 * randn(1, N);
    coefficients = [1 -3 -4];
    x1 = filter(coefficients, 1, w1);
    x1 = [x1, zeros(1, N2 - N)];
    X1 = fft(x1);
    Pxx = (1/N2) * abs(X1).^2;
    var = var + (abs(Pxx - Pxx_mean)).^2;
end
var = var / Mc;

subplot(2,2,3);
plot(freq(1:N2/2), var(1:N2/2), 'LineWidth', 2);
title('Variance - X1');
xlabel('Frequency (w)');
ylabel('var');
grid on;
xlim([0, pi]);

% MSE
MSE = var + B.^2;
subplot(2,2,4);
plot(freq(1:N2/2), MSE(1:N2/2), 'LineWidth', 2);
title('MSE - X1');
xlabel('Frequency (w)');
ylabel('MSE');
grid on;
xlim([0, pi]);

% Summary statistics (averages)
B_sum = sum(B(1:M)) / M;
var_sum = sum(var(1:M)) / M;
MSE_sum = sum(MSE(1:M)) / M;

disp(['B_sum = ' num2str(B_sum)]);
disp(['var_sum = ' num2str(var_sum)]);
disp(['MSE_sum = ' num2str(MSE_sum)]);
