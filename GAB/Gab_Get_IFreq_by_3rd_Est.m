function IFreq = Gab_Get_IFreq_by_3rd_Est(x, N_FFT, L, CR, method)
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
%   related paper: Abratkiewicz, K., & Gambrych, J. (2022). 
%   Real-time variants of vertical synchrosqueezing: Application 
%   to radar remote sensing. 
%   IEEE Journal of Selected Topics in Applied Earth Observations 
%   and Remote Sensing, 15, 1760-1774.

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
SdW = zeros(N_FFT, x.N);
Sd2W = zeros(N_FFT, x.N);
Sd3W = zeros(N_FFT, x.N);
StW = zeros(N_FFT, x.N);
St2W = zeros(N_FFT, x.N);
St2dW = zeros(N_FFT, x.N);
St2d2W = zeros(N_FFT, x.N);
StdW = zeros(N_FFT, x.N);
Std2W = zeros(N_FFT, x.N);

K = N_FFT;
mm = 1:N_FFT;
mm = mm(:);
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;

if strcmp(method, 'ptByPt')
    for n = 1:x.N
        k_min = min(n-1, round(K/2));
        k_max = min(x.N-n, round(K/2));
        k = ((-k_min):k_max).';
        w     = Gab_Gaussian_Window(k, L, 0, 0, CR);
        dw    = Gab_Gaussian_Window(k, L, 1, 0, CR);
        d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR);
        d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR);
        tw    = Gab_Gaussian_Window(k, L, 0, 1, CR);
        t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR);
        t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR);
        t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR);
        td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR);
        tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR);

        for m = 1:N_FFT
            S(m,n)      = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w     .* exp(A .* f_points(m) .* k));
            SdW(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* dw    .* exp(A .* f_points(m) .* k));
            Sd2W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* d2w   .* exp(A .* f_points(m) .* k));
            Sd3W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* d3w   .* exp(A .* f_points(m) .* k));
            StW(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw    .* exp(A .* f_points(m) .* k));
            StdW(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tdw   .* exp(A .* f_points(m) .* k));
            Std2W(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* td2w  .* exp(A .* f_points(m) .* k));
            St2W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2w   .* exp(A .* f_points(m) .* k));
            St2dW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2dw  .* exp(A .* f_points(m) .* k));
            St2d2W(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2d2w .* exp(A .* f_points(m) .* k));
        end

    end
elseif strcmp(method, 'FFT')

    k = -K/2:K/2-1;
    w     = Gab_Gaussian_Window(k, L, 0, 0, CR);
    dw    = Gab_Gaussian_Window(k, L, 1, 0, CR);
    d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR);
    d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR);
    tw    = Gab_Gaussian_Window(k, L, 0, 1, CR);
    t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR);
    t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR);
    t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR);
    td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR);
    tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR);
    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(0,K-1,N_FFT).';
    n_bins = 1:x.N;

    for n = 1:x.N
        tmp = signal(n_bins(n):n_bins(n)+K-1);
        S(:,n)      = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w, N_FFT));
        SdW(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw, N_FFT));
        Sd2W(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* d2w, N_FFT));
        Sd3W(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* d3w, N_FFT));
        StW(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw, N_FFT));
        StdW(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tdw, N_FFT));
        Std2W(:,n)  = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* td2w, N_FFT));
        St2W(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t2w, N_FFT));
        St2dW(:,n)  = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t2dw, N_FFT));
        St2d2W(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t2d2w, N_FFT));
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end
IFreq = mm - round(N_FFT.*(imag((Sd3W.*StdW.*St2W - StdW.*SdW.*St2d2W - Sd2W.*St2W.*Std2W - Sd3W.*St2dW.*StW + Sd2W.*St2d2W.*StW + SdW.*St2dW.*Std2W)./(Sd2W.*StdW.*St2W - Sd2W.*St2dW.*StW - SdW.*St2W.*Std2W + SdW.*St2d2W.*StW - StdW.*St2d2W.*S + St2dW.*Std2W.*S)))./2./pi);
end

