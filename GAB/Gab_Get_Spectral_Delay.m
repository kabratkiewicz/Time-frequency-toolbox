function SD = Gab_Get_Spectral_Delay(x, N_FFT, L, CR, method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)
%  CR - chirp rate of the window
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Płotka, M., Abratkiewicz, K., Malanowski, M., Samczyński, P., & Kulpa, K. 
%   (2020). The use of the reassignment technique in the time-frequency analysis 
%   applied in VHF-based passive forward scattering radar. Sensors, 20(12), 3434.

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
StW = zeros(N_FFT, x.N);

K = N_FFT;
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;

if strcmp(method, 'ptByPt')
    for n = 1:x.N
        k_min = min(n-1, round(K/2));
        k_max = min(x.N-n, round(K/2));
        k = ((-k_min):k_max).';
        w  = Gab_Gaussian_Window(k, L, 0, 0, CR);
        tw = Gab_Gaussian_Window(k, L, 0, 1, CR);
        for m = 1:N_FFT
            S(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w   .* exp(A .* f_points(m) .* k));
            StW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw  .* exp(A .* f_points(m) .* k));
        end
    end
elseif strcmp(method, 'FFT')
    k = -K/2:K/2-1;
    w  = Gab_Gaussian_Window(k, L, 0, 0, CR);
    tw = Gab_Gaussian_Window(k, L, 0, 1, CR);
    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(0,K-1,N_FFT).';
    n_bins = 1:x.N;
    for n = 1:x.N
        tmp = signal(n_bins(n):n_bins(n)+K-1);
        S(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w, N_FFT));
        StW(:,n)  = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw, N_FFT));
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end
nn = 1:x.N;
SD = nn + round(real(StW./S));
end
