function CR = Gab_Calculate_CR(x, N_FFT, L, CR_win, est, method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)
%  CR_win - chirp rate of the window
%  est - type of the estimator
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Samczyński, P., & Czarnecki, K. (2019). 
%   Radar signal parameters estimation using phase accelerogram in the 
%   time-frequency domain. IEEE Sensors Journal, 19(13), 5078-5085.

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
if ~exist('CR_win', 'var')
    CR_win = 0;
end
if ~exist('est', 'var')
    est = 3;
end
if ~exist('method', 'var')
    method = 'FFT';
end

K = N_FFT;
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;

if est == 1
    S   = zeros(N_FFT, x.N);
    SdW = zeros(N_FFT, x.N);
    StW = zeros(N_FFT, x.N);
    if strcmp(method, 'ptByPt')
        for n = 1:x.N
            k_min = min(n-1, round(K/2));
            k_max = min(x.N-n, round(K/2));
            k = ((-k_min):k_max).';
            w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
            dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
            tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
            for m = 1:N_FFT
                S(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w  .* exp(A .* f_points(m) .* k));
                SdW(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* dw .* exp(A .* f_points(m) .* k));
                StW(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw .* exp(A .* f_points(m) .* k));
            end
        end
    elseif strcmp(method, 'FFT')

        signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
        phase_shift = linspace(0,K-1,N_FFT).';
        n_bins = 1:x.N;
        k = -K/2:K/2-1;
        w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
        dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
        tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
        for n = 1:x.N
            tmp = signal(n_bins(n):n_bins(n)+K-1);
            S(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w,  N_FFT));
            SdW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw, N_FFT));
            StW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw, N_FFT));
        end
    else
        error('Wrong transform method. Available methods ptByPt or FFT');
    end

    nomK = real(SdW./S);
    denK = imag(StW./S);

    CR = (1/2/pi).*(nomK./denK);
elseif est == 2
    S   = zeros(N_FFT, x.N);
    SdW = zeros(N_FFT, x.N);
    Sd2W = zeros(N_FFT, x.N);
    SdtW = zeros(N_FFT, x.N);
    StW = zeros(N_FFT, x.N);
    if strcmp(method, 'ptByPt')
        for n = 1:x.N
            k_min = min(n-1, round(K/2));
            k_max = min(x.N-n, round(K/2));
            k = ((-k_min):k_max).';
            w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
            dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
            d2w = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
            tdw = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
            tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
            for m = 1:N_FFT
                S(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w   .* exp(A .* f_points(m) .* k));
                SdW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* dw  .* exp(A .* f_points(m) .* k));
                StW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw  .* exp(A .* f_points(m) .* k));
                Sd2W(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* d2w .* exp(A .* f_points(m) .* k));
                SdtW(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tdw .* exp(A .* f_points(m) .* k));
            end
        end
    elseif strcmp(method, 'FFT')
        K = round2even(K);
        k = -K/2:K/2-1;
        w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
        dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
        d2w = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
        tdw = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
        tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
        signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
        phase_shift = linspace(0,K-1,N_FFT).';
        n_bins = 1:x.N;
        for n = 1:x.N
            tmp =  signal(n_bins(n):n_bins(n)+K-1);
            S(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w,   N_FFT));
            SdW(:,n)  = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw,  N_FFT));
            StW(:,n)  = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw,  N_FFT));
            Sd2W(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* d2w, N_FFT));
            SdtW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tdw, N_FFT));
        end
    else
        error('Wrong transform method. Available methods ptByPt or FFT');
    end

    nomD = Sd2W.*S - (SdW).^2;
    denD = SdW.*StW - SdtW.*S;

    CR = (1/2/pi).*imag(nomD./denD);
elseif est == 3
    S   = zeros(N_FFT, x.N);
    SdW = zeros(N_FFT, x.N);
    St2W = zeros(N_FFT, x.N);
    SdtW = zeros(N_FFT, x.N);
    StW = zeros(N_FFT, x.N);
    if strcmp(method, 'ptByPt')
        for n = 1:x.N
            k_min = min(n-1, round(K/2));
            k_max = min(x.N-n, round(K/2));
            k = ((-k_min):k_max).';
            w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
            dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
            tdw = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
            tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
            t2w = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
            for m = 1:N_FFT
                S(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w   .* exp(A .* f_points(m) .* k));
                SdW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* dw  .* exp(A .* f_points(m) .* k));
                StW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw  .* exp(A .* f_points(m) .* k));
                St2W(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2w .* exp(A .* f_points(m) .* k));
                SdtW(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tdw .* exp(A .* f_points(m) .* k));
            end
        end
    elseif strcmp(method, 'FFT')
        K = round2even(K);
        k = -K/2:K/2-1;
        w  = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
        dw = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
        tdw = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
        tw = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
        t2w = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
        signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
        phase_shift = linspace(0,K-1,N_FFT).';
        n_bins = 1:x.N;

        for n = 1:x.N
            tmp = signal(n_bins(n):n_bins(n)+K-1);
            S(:,n)   = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w,  N_FFT));
            SdW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw, N_FFT));
            StW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw, N_FFT));
            St2W(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t2w, N_FFT));
            SdtW(:,n) = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tdw, N_FFT));
        end
    else
        error('Wrong transform method. Available methods ptByPt or FFT');
    end

    nomF = imag(SdtW./S - (SdW./S).*(StW./S));
    denF = real(St2W./S - (StW./S).^2);

    CR = -(1/2/pi).*(nomF./denF);
elseif est == 4
    S      = zeros(N_FFT, x.N);
    SdW    = zeros(N_FFT, x.N);
    Sd2W   = zeros(N_FFT, x.N);
    Sd3W   = zeros(N_FFT, x.N);
    StW    = zeros(N_FFT, x.N);
    St2W   = zeros(N_FFT, x.N);
    StdW   = zeros(N_FFT, x.N);
    Std2W  = zeros(N_FFT, x.N);
    St2dW  = zeros(N_FFT, x.N);
    St2d2W = zeros(N_FFT, x.N);
    if strcmp(method, 'ptByPt')
        for n = 1:x.N
            k_min = min(n-1, round(K/2));
            k_max = min(x.N-n, round(K/2));
            k = ((-k_min):k_max).';

            w     = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
            dw    = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
            d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
            d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR_win);
            tw    = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
            t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
            tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
            td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR_win);
            t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR_win);
            t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR_win);

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
        K = round2even(K);
        k = -K/2:K/2-1;
        w     = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
        dw    = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
        d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
        d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR_win);
        tw    = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
        t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
        tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
        td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR_win);
        t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR_win);
        t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR_win);
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

    CR = imag((-St2W.*Sd2W.^2 + St2dW.*Sd2W.*SdW + St2d2W.*S.*Sd2W - St2d2W.*SdW.^2 + Sd3W.*St2W.*SdW - Sd3W.*St2dW.*S)./(Sd2W.*StdW.*St2W - Sd2W.*St2dW.*StW - SdW.*St2W.*Std2W + SdW.*St2d2W.*StW - StdW.*St2d2W.*S + St2dW.*Std2W.*S))./2/pi;
   
else
    error('Wrong estimator type. Available estimators: 1 - K, 2 - D, 3 - F, 4 - in jerk');
end

