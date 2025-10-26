function [S] = my_stft(s, w_analysis, R, NFFT)
%MY_STFT Short-Time Fourier Transform
%   Computes the STFT of a signal s using the analysis window w_analysis,
%   hop size R, and FFT length NFFT.
%
% Inputs:
%   s           - input signal (column vector)
%   w_analysis  - analysis window (column vector)
%   R           - hop size between frames
%   NFFT        - FFT length (number of frequency bins)
%
% Output:
%   S           - STFT matrix (NFFT x Num_of_frames)

Win_len = length(w_analysis);  % Analysis window length
Num_of_frames = fix((length(s) - (Win_len - R)) / R);  % Number of frames
S = zeros(NFFT, Num_of_frames);  % Preallocate STFT matrix

for frame = 1:Num_of_frames
    % Extract and window the current segment
    segment = s((frame-1)*R + (1:Win_len)) .* w_analysis;
    
    % Compute FFT of the windowed segment
    temp_fft = fft(segment, Win_len);
    
    % If NFFT differs from Win_len, decimate/resize FFT to NFFT
    if NFFT ~= Win_len
        decimation_factor = ceil(Win_len / NFFT);
        temp_fft = downsample(temp_fft, decimation_factor);
        temp_fft = temp_fft(1:NFFT);  % Ensure length matches NFFT
    end
    
    % Apply phase shift for each frame
    for i = 1:NFFT
        S(i, frame) = temp_fft(i) * exp(-1j * 2 * pi / NFFT * (i-1) * (frame-1));
    end
end

end
