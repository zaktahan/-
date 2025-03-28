function s_hat_n = short_filter(S, R, L)
    % S - STFT of the signal
    % R - overlap samples (hop)
    % L - number of times the filter is being chopped
    
    [M, N] = size(S); 
    s_hat = zeros(M, N);
    Lw = L * R;
    s = 0:(Lw - 1);
    w_t = exp(-s / Lw); % as question defines w[s]
    
    if length(w_t) / L < length(S(1, :))
        w_t = [w_t, zeros(1, length(S(1, :)) - length(w_t))];
    end
 
    g = reshape(w_t, [], L)';
    [n, m] = size(g);
    G = zeros(n, m);
    
    % Constructing G
    for i = 1:n
        G(i, :) = fft(g(i, :));
    end
 
    Gn = zeros(size(s_hat));
    Gn(1:size(G, 1), 1:size(G, 2)) = G;
    s_hat_n = zeros(size(Gn));
    s_n = zeros(size(Gn));
    s_hat_n(1:size(s_hat, 1), 1:size(s_hat, 2)) = s_hat;
    s_n(1:size(S, 1), 1:size(S, 2)) = S;
 
    % Applying the filter
    for k = 1:M
        for n = 1:N
            temp = 0;
            for l = 1:L
                p = mod(l - 1, length(S(k, :))) + 1;
                temp = temp + Gn(l, mod(k - 1, length(Gn(k, :))) + 1) * s_n(k, p);
            end
            s_hat_n(k, n) = temp;
        end
    end
end

%% Time-varying filter
clear all
close all
clearvars

M = 512;
Lh = M;
Lg = M/2;
R = M/2;
Lf = M/2;
w_analysis = boxcar(Lh);
w_synthesis = boxcar(Lf);

s = cos(2*pi/512*97*(00:99999));fs=8000;

S = my__stft(s, w_analysis, R, M);
% L = size(S, 2); 
% 
% % Generate a time-varying frequency-domain filter
% Lw = L * Lg;
% s_values = 0:(Lw-1);
% w = exp(-s_values / Lw);  % Smooth exponential decay
% 
% % Reshape to match the STFT dimensions
% g_l_mat = reshape(w, Lg, L) / M;
% G_k_l = fft(g_l_mat, M);  % Get frequency response
% 
% % Apply filtering in the frequency domain (element-wise multiplication)
% S_filter = G_k_l .* S;
% 
% % Reconstruct the signal using ISTFT
% s_hat_filtered = real(my__istft(S_filter, w_synthesis, R));
% 
% % Normalize the output to prevent scaling issues
% s_hat_filtered = s_hat_filtered / max(abs(s_hat_filtered));
% 
% % Direct time-domain filtering for comparison
% s_filtered = filter(w, 1, s);
% s_filtered = s_filtered / max(abs(s_filtered));
s_hat_filtered = short_filter(S, R, 1)
% Play the reconstructed filtered signal
soundsc(s_hat_filtered, fs);

% Plot the signals
figure;

subplot(2, 1, 1);
plot(s_hat_filtered);
title('Filtered Signal using ISTFT (s\_hat\_filtered)');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(s_filtered);
title('Filtered Signal using FIR Filter (s\_filtered)');
xlabel('Time');
ylabel('Amplitude');
