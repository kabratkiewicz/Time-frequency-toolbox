function S = Gab_ASTFT(x, N_FFT, L,CR, method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - vector of a window duration parameter:  w0 * T, (default: 10)
%  CR - window chirp rate - scalar
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., & Samczyński, P. (2020, April). 
%   An adaptive spectrogram and accelerogram algorithm for electronic 
%   warfare applications. In 2020 IEEE International Radar Conference (RADAR) 
%   (pp. 530-535). IEEE.


signal = x.signal(:);

if ~exist('x', 'var')
    error('signal structure not defined');
end
if ~exist('N_FFT', 'var')
    N_FFT = 1024;
end
if ~exist('L', 'var')
    L = 10;
end
if ~exist('CR', 'var')
    CR = 0;
end
if ~exist('method', 'var')
    method = 'FFT';
end

S = zeros(N_FFT, x.N);
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;
K = N_FFT;
if strcmp(method, 'ptByPt')
    for n = 1:x.N
        k_min = min(n-1, round(K/2));
        k_max = min(x.N-n, round(K/2));
        k = ((-k_min):k_max).';
        w = Gab_Gaussian_Window(k, L(n), 0, 0, CR);
        for m = 1:N_FFT
            S(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w .* exp(A .* f_points(m) .* k));
        end
    end
elseif strcmp(method, 'FFT')

    k = -K/2:K/2-1;

    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(0,K-1,N_FFT).';
    n_bins = 1:x.N;
    for n = 1:x.N
        w = Gab_Gaussian_Window(k, L(n), 0, 0, CR);
        tmp = signal(n_bins(n):n_bins(n)+K-1);
        S(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w, N_FFT));
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end


end