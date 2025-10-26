% question1.m
% Spectrum plots for x1[n] and x2[n]
close all; clear; clc;

% Define the frequency range
w = linspace(0, pi, 2049);

% Theoretical PSD for x1[n]
sxx1 = (26 - 8*cos(2*w) + 18*cos(w)) / 26;

% Plot the Spectrum of x_1[n]
figure;
plot(w, sxx1);
title('Spectrum of x_1[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
grid on;

% Theoretical PSD for x2[n]
sxx2 = 0.51 ./ (1.49 - 1.4 * cos(w));
sxx2 = abs(sxx2);

% Plot the Spectrum of x_2[n]
figure;
plot(w, sxx2);
title('Spectrum of x_2[n]');
xlabel('Frequency (radians/sample)');
ylabel('Power/Frequency');
xlim([0 pi]);
grid on;
