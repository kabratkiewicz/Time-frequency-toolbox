function [S_mean, S_array] = Gab_Multitaper_STFT(x, N_FFT, L, M, method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)
%  M - number of tapers
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., & Samczyński, P. (2021, September). 
%   Multitaper time-frequency reassigned spectrogram in micro-Doppler radar 
%   signal analysis. In 2021 Signal Processing Symposium (SPSympo) (pp. 1-5). IEEE.

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
if ~exist('method', 'var')
    method = 'FFT';
end
K = N_FFT;
S_mean = zeros(N_FFT, x.N);
S_array = zeros(N_FFT, x.N, M);

f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
f_points = f_points(:);
A = -1j * 2*pi / N_FFT;

if strcmp(method, 'ptByPt')
    K = N_FFT;
    for i = 1:M
        for n = 1:x.N
            k_min = min(n-1, round(K/2));
            k_max = min(x.N-n, round(K/2));
            k = ((-k_min):k_max).';
            [w,~,~] = hermf(K,M,L) ; 
            for m = 1:N_FFT
                S_mean(m,n, i) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w(:,m).* exp(A .* f_points(m) .* k));
            end
        end
    end
elseif strcmp(method, 'FFT')
   
    [w,~,~] = hermf(K,M,L) ; 
    w = w.';
    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(0,K-1,N_FFT).';
    n_bins = 1:x.N;
    for m = 1:M
        for n = 1:x.N
             tmp = signal(n_bins(n):n_bins(n)+K-1);
             S_array(:,n,m) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w(:,m), N_FFT));
        end
    end
    S_mean = mean(abs(S_array).^2,3);
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end

end

