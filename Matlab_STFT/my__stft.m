%stft is a short time fourier transformation-we take a small portion of the
%whole signal dn fft it:
function [S] = my_stft(s, w_analysis, R, NFFT)
%parameters:

    Win_len = length(w_analysis);% the winddow we multiply the part of the signal befor fft'ing
    Num_of_frames = fix((length(s) - (Win_len - R)) / R);% num of the parts we will work with
    S = zeros(NFFT, Num_of_frames);

for frame = 1:Num_of_frames
    % Initialize x to store the windowed segment
    x = zeros(Win_len, 1);
    
    % multiplying the part of the signal with the window:
    for k = 1:Win_len
        x(k) = w_analysis(k) * s((frame - 1) * R + k);
    end
    
    % fft'ing and inserting to our output S:
    temp = fft(x, Win_len);
        for i = 1:NFFT
            decimated = downsample(temp, ceil(Win_len / NFFT));
            S(i, frame) = decimated(i) * exp(-1j * 2 * pi / NFFT * (i-1) * (frame-1));%shift
        end
end
