function S = Gab_CT(x, N_FFT, L, gamma_K, CR, method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)
%  CR - vector of the window chirp rate 
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K. (2020), Double-adaptive chirplet transform for 
%   radar signature extraction. IET Radar Sonar Navig., 14: 1463-1474. 
%   https://doi.org/10.1049/iet-rsn.2020.0084

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
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end
if ~exist('CR', 'var')
 CR = zeros(1,x.N);
end
if ~exist('method', 'var')
 method = 'FFT';
end

S = zeros(N_FFT, x.N);
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;

if strcmp(method, 'ptByPt')
    for n = 1:x.N
        K = 2 * L(n) * sqrt(2*log(1/gamma_K));  %% window length in samples
        k_min = min(n-1, round(K/2));
        k_max = min(x.N-n, round(K/2));
        k = ((-k_min):k_max).';
        w = Gab_Gaussian_Window(k, L, 0, 0, CR(n));
        for m = 1:N_FFT
            S(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w .* exp(A .* f_points(m) .* k));
        end
    end
elseif strcmp(method, 'FFT') 
 K = N_FFT;
    k = -K/2:K/2-1;
     
    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(0,K-1,N_FFT).';
    n_bins = 1:x.N;
    
    for n = 1:x.N
        tmp = signal(n_bins(n):n_bins(n)+K-1);
        w = Gab_Gaussian_Window(k, L, 0, 0, CR(n));
        S(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w, N_FFT));
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end


end


