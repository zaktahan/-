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

[s,fs] = audioread('fdmy0-sx297.wav');
%s = cos(2*pi/512*97*(00:99999));

S = my_stft(s, w_analysis, R ,M);
L = size(S, 2); 
s_hat = real(my_istft(S, w_synthesis, R));

Lw = L*Lg;
s_values = 0:(Lw-1);
w = exp(-s_values/Lw);

g_l_mat = reshape(w,Lg, L)/M;
G_k_l = fft(g_l_mat, M);

S_filter = zeros(size(S, 1),size(S, 2));
for k=1:M
    S_filter(k,:)= filter(G_k_l(k,:), 1, S(k,:));
end 

s_hat_filtered =  real(my_istft(S_filter, w_synthesis, R));

s_filtered = filter(w,1,s);
soundsc(s_hat_filtered)
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