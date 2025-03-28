clear all
close all
% parameters
n = 100;
m = 5;
X = rand(n, m);
beta = [2; -1; 3; 0.5; 4];
sigma = 0.1;
epsilon = sigma * randn(n, 1);
% define Y
Y = X * beta + epsilon;

%beta^ = (X' * X)^-1 * X' * Y
beta_hat = (X' * X) \ (X' * Y);

% compare b,b^ as function of sigmas
disp('Original beta values:');
disp(beta);
disp('Estimated beta_hat values:');
disp(beta_hat);
sigmas = [0.1, 0.5, 1, 2];
errors = zeros(length(sigmas), 1);
%for each sigma
for i = 1:length(sigmas)
    epsilon = sigmas(i) * randn(n, 1); 
    Y = X * beta + epsilon; 
    beta_hat = (X' * X) \ (X' * Y); 
    errors(i) = norm(beta - beta_hat); 
end

% plot
figure;
plot(sigmas, errors, '-o');
xlabel('\sigma (Standard Deviation of Noise)');
ylabel('Error between beta and beta_hat');
title('Effect of Noise on Estimation Accuracy');
grid on;
