function [S, iFreq, gDel] = Gab_Multitaper_Reass_Vectors(x, N_FFT, L, M)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)
%  M - number of tapers
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

S = zeros(N_FFT, x.N, M);
iFreq = zeros(N_FFT, x.N, M);
gDel = zeros(N_FFT, x.N, M);
Sdw = zeros(N_FFT, x.N);
Stw = zeros(N_FFT, x.N);

f_points = 1:N_FFT;
f_points = f_points(:);

K = N_FFT;
[w,dw,~] = hermf(K,M,L) ;
w = w.';
dw = dw.';
tw = linspace(-K/2, K/2, K).'.*w;
signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
phase_shift = linspace(0,K-1,N_FFT).';
n_bins = 1:x.N;
for m = 1:M
    for n = 1:x.N
        tmp = signal(n_bins(n):n_bins(n)+K-1);
        S(:,n,m)    = exp(1i * pi * phase_shift) .*fftshift(fft(tmp .* w(:,m), N_FFT));
        Stw(:,n)  = exp(1i * pi * phase_shift) .*fftshift(fft(tmp .* tw(:,m), N_FFT));
        Sdw(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw(:,m), N_FFT));
    end
    iFreq(:,:,m) = f_points - round(imag(N_FFT .* (Sdw./S(:,:,m))/2/pi));
    gDel(:,:,m) = n_bins + round(real( (Stw./S(:,:,m))));
end

end

