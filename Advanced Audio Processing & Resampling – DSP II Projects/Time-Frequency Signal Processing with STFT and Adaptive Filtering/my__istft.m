function [s_hat] = my_istft(S, w_synthesis, R)
%MY_ISTFT Inverse Short-Time Fourier Transform (ISTFT)
%   Reconstructs time-domain signal from STFT matrix S
%   using the provided synthesis window w_synthesis and hop size R.

% Inputs:
%   S           - STFT matrix (frequency bins x frames)
%   w_synthesis - synthesis window (column vector)
%   R           - hop size between frames
%
% Output:
%   s_hat       - reconstructed time-domain signal

[Num_of_frequency_bands, Num_of_frames] = size(S);
Win_len = length(w_synthesis);                     % Synthesis window length

% Preallocate output signal
s_hat = zeros(R * Num_of_frames + Win_len - R, 1);

% Loop over each frame
for k = 1:Num_of_frames
    current_freq = S(:, k);
    
    % Apply phase shift for current frame
    for n = 1:Num_of_frequency_bands
        current_freq(n) = current_freq(n) * exp(1j * 2 * pi / Num_of_frequency_bands * (n-1) * (k-1));
    end
    
    % Compute inverse FFT
    temp = ifft(current_freq);
    
    % Apply synthesis window with wrapping if necessary
    s_reconstructed = zeros(Win_len, 1);
    for n = 1:Win_len
        idx = mod(n-1, Num_of_frequency_bands) + 1;
        s_reconstructed(n) = temp(idx) * w_synthesis(n);
    end
    
    % Overlap-add to the output signal
    start_index = (k - 1) * R + 1;
    s_hat(start_index:start_index + Win_len - 1) = s_hat(start_index:start_index + Win_len - 1) + s_reconstructed;
end

end
