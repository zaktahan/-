%% -------------------- Parameters Initialization --------------------
Fs = 8000;       % Original sampling frequency
L  = 4;          % Upsampling factor
M  = 3;          % Downsampling factor

% Filter order: depends on max(L,M) to ensure a narrow and steep filter
N = 20 * max(L, M);  % Filter order

% Low-pass filter (LPF) for anti-aliasing and interpolation
if L ~= 1 || M ~= 1
    h = L * fir1(N-1, 1/max(L, M));
end

%% -------------------- Load Audio Files --------------------
% NOTE: Change the file paths below to match the location on your computer
[doors, Fs_orig0]     = audioread('C:\Users\zakta\OneDrive\שולחן העבודה\doors.wav');
[noAnswer, Fs_orig1]  = audioread('C:\Users\zakta\OneDrive\שולחן העבודה\fdaw0-si1406.wav');
[rockNroll, Fs_orig2] = audioread('C:\Users\zakta\OneDrive\שולחן העבודה\fdmy0-sx297.wav');
[mugs, Fs_orig3]      = audioread('C:\Users\zakta\OneDrive\שולחן העבודה\medr0-si744.wav');
[laws, Fs_orig4]      = audioread('C:\Users\zakta\OneDrive\שולחן העבודה\mewm0-si718.wav');

% Choose which audio file to process
signal = doors;

%% -------------------- Resampling --------------------
resampled = resample_signal(signal, L, M, N, Fs, h);

% Normalization options
resampled_signal_Lnorm = resampled;          % Default normalization by L
resampled_signal_Mnorm = (M/L) * resampled;  % Normalization by M to match original magnitude

%% -------------------- FFT and Plotting --------------------
X         = fft(signal, length(signal));
Y_Mnorm   = fft(resampled_signal_Mnorm, length(resampled));
Y_Lnorm   = fft(resampled_signal_Lnorm, length(resampled));

f_axis_x = (0:length(X)-1) * (Fs / length(X));
f_axis_y = (0:length(Y_Mnorm)-1) * (Fs / length(Y_Mnorm));

figure;
subplot(4,1,1)
plot(f_axis_x, abs(X));
title('Spectrum of Original Signal x[n]');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]);

subplot(4,1,2)
plot(f_axis_y, abs(Y_Mnorm));
title('Spectrum of Output Signal y[n] Mnorm');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]);

subplot(4,1,3)
plot(f_axis_y, abs(Y_Lnorm));
title('Spectrum of Output Signal y[n] Lnorm');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]);

% Comparison with MATLAB built-in resample
if L == 1 && M == 1
    y_builtin = resample(signal, L, M);
else
    y_builtin = resample(signal, L, M, h);
end

Y_builtin       = fft(y_builtin, length(y_builtin));
f_axis_y_builtin = (0:length(Y_builtin)-1) * (Fs / length(Y_builtin));

subplot(4,1,4)
plot(f_axis_y_builtin, abs(Y_builtin));
title('Spectrum of Output From MATLAB Resample');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 Fs/2]);

%% -------------------- Audio Playback --------------------
soundsc(signal, Fs);          pause(4);
soundsc(resampled_signal_Lnorm, Fs*L/M); pause(4);
soundsc(resampled_signal_Lnorm, Fs*L/M); pause(4);
soundsc(y_builtin, Fs*L/M);

%% -------------------- Resample Function --------------------
function resampled_signal = resample_signal(signal, L, M, N, Fs, h)
    % Handle trivial case where no resampling is needed
    if L == 1 && M == 1
        resampled_signal = signal;
        return;
    end

    signal_length = length(signal);
    output_len    = floor(signal_length * L / M);  
    Q             = ceil(N / L);  % Polyphase buffer length

    n_prime = 0:output_len-1;      % Output sample indices
    I       = (0:(Q-1))';          % Polyphase branch indices

    % Indices matrices for signal and filter
    buf_idx    = floor(n_prime * M / L);
    signal_idx = buf_idx - I;
    filter_idx = I * L + mod(n_prime * M, L);
    filter_idx = mod(filter_idx, length(h)); % Wrap-around filter indices

    % Replace invalid indices
    valid = (signal_idx >= 1) & (signal_idx <= length(signal));
    signal_idx(~valid) = 1;

    % Build signal and filter matrices
    signal_values = signal(signal_idx);
    filter_values = h(filter_idx + 1);

    % Polyphase multiplication and sum
    polyphase_matrix      = signal_values .* filter_values;
    polyphase_matrix(~valid) = 0;

    resampled_signal = sum(polyphase_matrix, 1);
end
