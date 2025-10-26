%% Linear Regression Parameter Estimation
clear all;
close all;

%% Parameters
n = 100;
m = 5;
X = rand(n, m);
beta = [2; -1; 3; 0.5; 4];
sigma = 0.1;
epsilon = sigma * randn(n, 1);

%% Generate Target Variable
Y = X * beta + epsilon;

%% Estimate Coefficients (β̂ = (X'X)⁻¹ X'Y)
beta_hat = (X' * X) \ (X' * Y);

%% Display Results
disp('Original beta values:');
disp(beta);
disp('Estimated beta_hat values:');
disp(beta_hat);

%% Compare β and β̂ for Different Noise Levels
sigmas = [0.1, 0.5, 1, 2];
errors = zeros(length(sigmas), 1);

for i = 1:length(sigmas)
    epsilon = sigmas(i) * randn(n, 1);
    Y = X * beta + epsilon;
    beta_hat = (X' * X) \ (X' * Y);
    errors(i) = norm(beta - beta_hat);
end

%% Plot Error vs Noise Level
figure;
plot(sigmas, errors, '-o', 'LineWidth', 1.5);
xlabel('\sigma (Standard Deviation of Noise)');
ylabel('Error between \beta and \betâ');
title('Effect of Noise on Estimation Accuracy');
grid on;
