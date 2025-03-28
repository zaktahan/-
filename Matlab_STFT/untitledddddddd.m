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
%s = cos(2*pi/512*97*(0:99999));

% Compute STFT
S = my__stft(s, w_analysis, R ,M);
L = 1%size(S, 2); 

% Compute the filter in frequency domain
Lw = L * Lg;
s_values = 0:(Lw-1);
w = exp(-s_values / Lw);
%%
g_l_mat = reshape(w, Lg, L) / M;
figure;
hold on;
for i = 1:size(g_l_mat, 2) % Iterate over columns (5 in total)
    plot(s_values((i-1)*Lg + (1:Lg)), g_l_mat(:, i)); % Use correct indexing
end
hold off;
legend("Column 1", "Column 2", "Column 3", "Column 4", "Column 5");
xlabel("Index");
ylabel("Amplitude");
title("Plot of Each Column of g_l_mat");
grid on;

%%
G = zeros(M,L);
G = fft(g_l_mat,M)/M;
[k, n] = size(S); % Extract the dimensions of S

% Calculate the new number of time samples
new_n = n + R * (L - 1);

% Zero-pad the matrix S to the new size
S_padded = [S, zeros(k, new_n - n)]; % Pad with zeros along the second dimension
S=S_padded;
S_filtered = zeros(M,size(S,2));
S_delayed=zeros(M,size(S,2));
for k=1:M%l=1:L
    for l=1:L%k=1:M
        S_delayed=delay(S,(l-1)*R);
        S_filtered(k,:) = S_filtered(k,:)+G(k,l)*S_delayed(k,:);
    end
end
%%
x=real(my__istft(S_filtered, w_synthesis, R));
if L > 1
    % Find the first non-zero index
    start_idx = find(x ~= 0, 1, 'first'); % First non-zero index

    % Find the next zero index after the first signal
    end_idx = find(x(start_idx:end) == 0, 1, 'first'); 

    if isempty(end_idx)
        % If there is no zero after the first signal, take the entire remaining part
        first_signal = x(start_idx:end);
        signal_start = length(x) + 1; % No second signal
    else
        % Extract the first continuous signal part
        first_signal = x(start_idx:start_idx + end_idx - 2); 
        signal_start = start_idx + end_idx; % Start after the first signal
    end

    % Find the second non-zero segment
    second_start_idx = find(x(signal_start:end) ~= 0, 1, 'first'); 

    if ~isempty(second_start_idx)
        second_start_idx = signal_start + second_start_idx - 1;  % Adjust the index to match the full vector
        second_end_idx = find(x(second_start_idx:end) == 0, 1, 'first');

        if isempty(second_end_idx)
            % If there is no zero after the second signal, take the entire remaining part
            second_signal = x(second_start_idx:end);
            signal_start = length(x) + 1; % No third signal
        else
            % Extract the second continuous signal part
            second_signal = x(second_start_idx:second_start_idx + second_end_idx - 2);
            signal_start = second_start_idx + second_end_idx; % Start after the second signal
        end
    else
        second_signal = [];  % If there's no second signal, return empty
        signal_start = length(x) + 1; % No third signal
    end

    x = second_signal;
end


%%
% Compare with time-domain filtering
s_filtered = filter(w, 1, s);

% Play filtered signal
soundsc(s_filtered, fs);

% Plot the signals
figure;

subplot(3, 1, 1);
plot(x);
title('Filtered Signal using ISTFT (s\_hat\_filtered)');
xlabel('Time');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(s_filtered);
title('Filtered Signal using FIR Filter (s\_filtered)');
xlabel('Time');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(s);
title('Filtered Signal using FIR Filter (s)');
xlabel('Time');
ylabel('Amplitude');