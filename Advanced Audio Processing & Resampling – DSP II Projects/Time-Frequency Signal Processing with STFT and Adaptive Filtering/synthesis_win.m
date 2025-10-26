function ws = synthesis_win(w, L)
% SYNTHESIS_WIN Computes the synthesis window for WOLA (Weighted Overlap-Add)
%
% ws = synthesis_win(w, L)
%
% Inputs:
%   w - analysis window vector (Nx1)
%   L - hop size (jump between windows)
%
% Output:
%   ws - synthesis window (Nx1)

N = length(w);
ws = w;                 % initialize synthesis window
w2 = zeros(N, 1);       % temporary vector to accumulate squared overlaps

for n = 0:N-1
    qn = floor(n/L);         % max number of backward overlaps
    qm = floor((N-1-n)/L);   % max number of forward overlaps
    for q = -qn:qm
        w2(n+1) = w2(n+1) + w(q*L + n + 1)^2;  % sum squared overlapping windows
    end
end

% compute synthesis window
ws = flipud(w ./ w2);

end
