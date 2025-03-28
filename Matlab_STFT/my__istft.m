function [s_hat] = my_istft(S, w_synthesis, R)

% Parameters:
[Num_of_frequency_bands, Num_of_frames] = size(S);
Win_len = length(w_synthesis);% Length of synthesis window-we will multiply the ifft result with it in the output
s_hat = zeros(R * Num_of_frames + Win_len - R, 1); % lenght of the resulted signal
% istft:
for k = 1:Num_of_frames% for every frequency band:
    s_reconstructed = zeros(Win_len, 1);
    current_freq = S(:, k);
    for n = 1:Num_of_frequency_bands
        current_freq(n) = current_freq(n) * exp(1j * 2 * pi / Num_of_frequency_bands * (n-1) * (k-1));%shift
    end
    % Compute inverse FFT for the current frame
    temp = ifft(current_freq);
    %multiplying the result with the synthesis window:the mod function is
    %becouse of the wrapping
    for n = 1:Win_len
        s_reconstructed(n) = temp(mod(n-1, Num_of_frequency_bands) + 1) * w_synthesis(n);
    end
%each time we will take the resulted band to the output-in jumps of R
     start_index = (k - 1) * R + 1;
    for i = 1:Win_len
        s_hat(start_index + i - 1) = s_hat(start_index + i - 1) + s_reconstructed(i);
    end
end       
end
