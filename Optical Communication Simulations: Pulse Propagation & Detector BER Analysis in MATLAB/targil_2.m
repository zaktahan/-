%% Optical Communication Project - Exercise 2

clear; clc; close all;

%% =======================
%% Q1: Received power for a given BER
%% =======================
delta_f = 1e10;   % 10 Gbps bandwidth [Hz]
T       = 300;    % [K] temperature
R_L     = 50;     % load resistance [Ohm]
BER     = 1e-9;   % target BER
R       = 1;      % [W/A] responsivity
Fn      = 3;      % noise figure

K_B = 1.38e-23;   % Boltzmann constant
q   = sqrt(2) * erfcinv(2 * BER); % BER = 1/2 * erfc(q/sqrt(2))
p_rec = (q/R) * sqrt(4*K_B*T*Fn*delta_f/R_L); % received power [W]

fprintf('Q1: Received power at BER=10^-9 = %.4g W\n', p_rec);

%% =======================
%% Q2: Plot received power vs BER
%% =======================
q_vals = linspace(2, 8, 1e3);
BER_vals = 0.5 * erfc(q_vals / sqrt(2));
p_rec_vals = q_vals .* sqrt(4*K_B*T*Fn*delta_f/R_L) / R;

figure;
plot(p_rec_vals, BER_vals, 'LineWidth', 1.5);
grid on;
xlabel('p_{rec} [W]');
ylabel('BER');
title('Received Power vs BER');

%% =======================
%% Q3: Simulation of PIN Detector
%% =======================
BER_target = 1e-3;
[~, idx] = min(abs(BER_vals - BER_target));
P_rec3 = p_rec_vals(idx);
P_in3  = 2 * P_rec3;

fprintf('Q3: Input power at BER=10^-3 = %.4g W\n', P_in3);

N_bits = 1e5;
symbols = randi([0, 1], 1, N_bits);

I_3 = symbols .* R * 2 * P_rec3;       % ideal current
I_d3 = P_in3 / 2;                       % decision threshold
sig_3 = sqrt((4*K_B*T/R_L)*Fn*delta_f);% thermal noise std dev

I_3_noise = I_3 + normrnd(0, sig_3, 1, N_bits);

% Count errors
errors = sum((symbols==0 & I_3_noise>I_d3) | (symbols==1 & I_3_noise<I_d3));
fprintf('Q3: Number of errors (PIN detector) = %d\n', errors);

% Plot histogram
figure;
histogram(I_3_noise, 100);
xlabel('I_{noise} [A]');
ylabel('Number of events');
title('Histogram of PIN Detector Current');
xline(I_d3, '--', 'I_d');

%% =======================
%% Q4 & Q5: Avalanche Photodiode Detector
%% =======================
M   = 10;       % multiplication factor
k_a = 0.7;      % excess noise factor parameter
F_a = k_a*M + (1-k_a)*(2-1/M);
e   = 1.6e-19;  % electron charge [C]

% Mean currents
I_0 = 0;
I_1 = 2*M*R*P_rec3;

% Noise std dev
sig_0 = sqrt((4*K_B*T/R_L)*Fn*delta_f);
sig_1 = sqrt((4*K_B*T/R_L)*Fn*delta_f + 4*e*M^2*R*P_rec3*F_a*delta_f);

% Decision current
I_d4 = (sig_0*I_1 + sig_1*I_0) / (sig_0 + sig_1);

% Simulate detector output with noise
I_4_noise = zeros(1, N_bits);
symbols   = randi([0, 1], 1, N_bits);
errors2   = 0;

for i = 1:N_bits
    if symbols(i) == 0
        I_4_noise(i) = I_0 + normrnd(0, sig_0);
        if I_4_noise(i) > I_d4
            errors2 = errors2 + 1;
        end
    else
        I_4_noise(i) = I_1 + normrnd(0, sig_1);
        if I_4_noise(i) < I_d4
            errors2 = errors2 + 1;
        end
    end
end

fprintf('Q4/5: Number of errors (APD detector) = %d\n', errors2);

% Plot histogram
figure;
histogram(I_4_noise, 100);
xlabel('I_{noise} [A]');
ylabel('Number of events');
title('Histogram of APD Detector Current');
xline(I_d4, '--', 'I_d');
