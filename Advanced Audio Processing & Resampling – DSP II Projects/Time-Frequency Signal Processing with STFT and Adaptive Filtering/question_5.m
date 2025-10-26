%% Time-Varying Filter - STFT-based Filtering
clc; clear; close all;

%% Parameters
M = 512;         % Number of frequency bins
Lh = M;          % Analysis window length
Lg = M/2;        % Filter length
R = M/2;         % Hop size
Lf = M/2;        % Synthesis window length
w_analysis = boxcar(Lh);   % Analysis window
w_synthesis = boxcar(Lf);  % Synthesis window

%% Load Signal
[s, fs] = audioread('fdmy0-sx297.wav');
% s = cos(2*pi/512*97*(0:99999)); % Optional test signal

%% Compute STFT
S = my__stft(s, w_analysis, R, M);
L = 1; % Number of filter segments (can be adjusted)

%% Construct Time-Varying Filter
Lw = L * Lg;
s_values = 0:(Lw-1);
w = exp(-s_values / Lw);

% Reshape filter into matrix g_l_mat
g_l_mat = reshape(w, Lg, L) / M;

% Plot columns of g_l_mat
figure;
hold on;
for i = 1:size(g_l_mat, 2)
    idx_range = (i-1)*Lg + (1:Lg);
    plot(s_values(idx_range), g_l_mat(:, i));
end
hold off;
legend(arrayfun(@(x) sprintf('Column %d', x), 1:size(g_l_mat,2), 'UniformOutput', false));
xlabel('Index');
ylabel('Amplitude');
title('Columns of g\_l\_mat');
grid on;

%% FFT of filter segments
G = fft(g_l_mat, M)/M;

%% Zero-pad STFT if needed
[k, n] = size(S);
new_n = n + R * (L - 1);
S = [S, zeros(k, new_n - n)];

%% Apply Time-Varying Filter
S_filtered = zeros(M, size(S,2));
for l = 1:L
    S_delayed = delay(S, (l-1)*R);
    for k = 1:M
        S_filtered(k,:) = S_filtered(k,:) + G(k,l) * S_delayed(k,:);
    end
end

%% Reconstruct Signal using ISTFT
x = real(my__istft(S_filtered, w_synthesis, R));

%% Extract second non-zero segment if L > 1
if L > 1
    start_idx = find(x ~= 0, 1, 'first');
    end_idx = find(x(start_idx:end) == 0, 1, 'first');
    if isempty(end_idx)
        first_signal = x(start_idx:end);
        signal_start = length(x)+1;
    else
        first_signal = x(start_idx:start_idx + end_idx - 2);
        signal_start = start_idx + end_idx;
    end

    second_start_idx = find(x(signal_start:end) ~= 0, 1, 'first');
    if ~isempty(second_start_idx)
        second_start_idx = signal_start + second_start_idx - 1;
        second_end_idx = find(x(second_start_idx:end) == 0, 1, 'first');
        if isempty(second_end_idx)
            second_signal = x(second_start_idx:end);
        else
            second_signal = x(second_start_idx:second_start_idx + second_end_idx - 2);
        end
    else
        second_signal = [];
    end
    x = second_signal;
end

%% Compare with Standard FIR Filtering
s_filtered = filter(w, 1, s);

%% Playback
soundsc(s_filtered, fs);

%% Plot Results
figure;
subplot(3,1,1);
plot(x);
title('Filtered Signal using ISTFT (S\_hat\_filtered)');
xlabel('Time');
ylabel('Amplitude');

subplot(3,1,2);
plot(s_filtered);
title('Filtered Signal using FIR Filter (s\_filtered)');
xlabel('Time');
ylabel('Amplitude');

subplot(3,1,3);
plot(s);
title('Original Signal (s)');
xlabel('Time');
ylabel('Amplitude');
