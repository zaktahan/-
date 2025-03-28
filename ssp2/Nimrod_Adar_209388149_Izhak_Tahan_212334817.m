% initializing parameters
Fs = 8000; % original sampling frequency
L = 4; % upsampling factor 
M = 3; % downsampling factor 

% filter order is dependent on max(L,M) since we would want the filter to be
% narrow and steep enough in frequency so we wont get leakage from aliased
% copies of the signal into the output, taking max(L,M) tackles both the
% case of high ratio in favor of M and the case of high ratio in favor of L
% the 20 factor has more to do with the original sampling frequency of 8000
% (very) generally speaking a filter of order 20 would suffice for speech filtering 
% and would be around the smaller values which retain enough quality which
% would be preferable for higher efficiancy.
N = 20*max(L,M); 

% LPF, which is dervied from the decimation's anti-aliasing filter and the interpolation filter
if L ~= 1 || M ~= 1
    h = L*fir1(N-1, 1/max(L, M));
end
%handling the case in which theres 1 at the cutoff and fir1 cant opperate 

% audio files
[doors, Fs_orig0] = audioread('/Users/nimo/Documents/DSP2_mat_Ex1 wav/doors.wav');
[noAnswer, Fs_orig1] = audioread('/Users/nimo/Documents/DSP2_mat_Ex1 wav/fdaw0-si1406.wav');
[rockNroll, Fs_orig2] = audioread('/Users/nimo/Documents/DSP2_mat_Ex1 wav/fdmy0-sx297.wav');
[mugs, Fs_orig3] = audioread('/Users/nimo/Documents/DSP2_mat_Ex1 wav/medr0-si744.wav');
[laws, Fs_orig4] = audioread('/Users/nimo/Documents/DSP2_mat_Ex1 wav/mewm0-si718.wav');

%choose audio file
signal = doors;
%apply the resampling function
resampled=resample_signal(signal, L, M, N, Fs, h);
    
    %defult normalisation by L is applied by the filter as shown in class
    %this gives similar magnitude to matlab's resample
    resampled_signal_Lnorm = resampled;

    %normlizing by M instead of by L get us the same relative magnitude 
    %as the original as expected from the theory (since only decimation effects the magnitude)
    resampled_signal_Mnorm = (M/L)*resampled;

    %plotting FFT of original, resampled_signal_Lnorm,
    %resampled_signal_Mnorm and matlab's built in resample for comparison.
    X = fft(signal, length(signal));
    Y_Mnorm = fft(resampled_signal_Mnorm, length(resampled));
    Y_Lnorm = fft(resampled_signal_Lnorm, length(resampled));

    f_axis_x = (0:length(X)-1) * (Fs / length(X));
    f_axis_y = (0:length(Y_Mnorm)-1) * (Fs / length(Y_Mnorm));

    figure;
    subplot(4, 1, 1);
    plot(f_axis_x, abs(X)); 
    title('Spectrum of Original Signal x[n]');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 Fs/2]);

    subplot(4, 1, 2);
    plot(f_axis_y, abs(Y_Mnorm)); 
    title('Spectrum of Output Signal y[n''] Mnorm');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 Fs/2]);

    subplot(4, 1, 3);
    plot(f_axis_y, abs(Y_Lnorm)); 
    title('Spectrum of Output Signal y[n''] Lnorm');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 Fs/2]);


%comparing between our function and matlab's built-in function
%using the same filter they give the same results
if L == 1 && M == 1
    y_builtin = resample(signal, L, M);
else
    y_builtin = resample(signal, L, M, h);
end
%handling the case in which theres 1 at the cutoff and fir1 cant opperate 

    Y_builtin = fft(y_builtin, length(y_builtin));
    f_axis_y_builtin = (0:length(Y_builtin)-1) * (Fs / length(Y_builtin));

    subplot(4, 1, 4);
    plot(f_axis_y_builtin, abs(Y_builtin));
    title('Spectrum of Output From Matlabs Resample');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    xlim([0 Fs/2]);

    % play back the original and resampled signals
    soundsc(signal, Fs);
    pause(4);
    soundsc(resampled_signal_Lnorm, Fs * L / M);
    pause(4);
    soundsc(resampled_signal_Lnorm, Fs * L / M);
    pause(4);
    soundsc(y_builtin, Fs * L / M);

 
function resampled_signal = resample_signal(signal, L, M, N, Fs, h)
    %optional - making sure the signal is a column vector
    %signal = signal(:);

    %optional - zero-padding if we dont want the original clicks in the end
    %padding_length = ceil(length(signal) * (L / M));
    %signal = [signal; zeros(padding_length, 1)];

%handling the case in which theres 1 at the cutoff and fir1 cant opperate 
%the function simply returns the signal unchanged
if L == 1 && M == 1
        resampled_signal = signal;
        return;
end

    %initialize parameters 
    signal_length = length(signal);
    output_len = floor(signal_length * L / M);  
    Q = ceil(N / L); % buffer length
   
%each output sample y[n'] is a sum over I from 0 to Q-1 of
%x(floor(n'M/L+I)*h((Q-1-I)*L+n'MmodL) so we would first build an
%indicies matrixes, one for the indicies of x and one for the indicies
%of h. (it was easier for us to write the polyphase directly through h's 
%indices and not go through the regular p notation).

n_prime = 0:output_len-1; % output sample indices, which would be the column indices
I = (0:(Q-1))'; % the branchs indices which would be row indices 

% signal idx matrix with values floor(n'M/L - I)
buf_idx = floor(n_prime * M / L); 
signal_idx = buf_idx - I; % input index matrix (Q x output_len)

% filter idx matrix with values (IL + (n'M mod L))
filter_idx = I * L + mod(n_prime * M, L); % filter index matrix (Q x output_len)
filter_idx = mod(filter_idx, length(h)); % wrap filter indices if necessary to stay within bounds

% replace invalid signal indices (negative or out of bounds) with 1 (arbitrary)
valid = (signal_idx >= 1) & (signal_idx <= length(signal));
signal_idx(~valid) = 1;  

% build matrices of signal values and filter coefficients using index matrices
signal_values = signal(signal_idx); 
filter_values = h(filter_idx + 1); 

% multipling the value matrices would produce a matrix with the values
%x(floor(n'M/L+I)*h((Q-1-I)*L+n'MmodL)

polyphase_matrix = signal_values .* filter_values;
polyphase_matrix(~valid) = 0; % ignore invalid contributions

%the output vector is obtained by summing all rows (the different I values)
%for each n', essentially summing the matrix's rows

resampled_signal = sum(polyphase_matrix, 1);

    
end


 