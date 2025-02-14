function [S, HSS, CR] = Gab_TF_H2_Synchrosqueezing(x, N_FFT, L, CR_win, IFE_method, CRE_method, Fourier_method)
%% x - structure of the signal and its params, fields:
%      .signal - vector of the signal samples
%      .fs     - sampling rate
%      .N      - signal length in samples
%      .T      - signal duration in [s]
%  N_FFT - number of frequency bins
%  L - window duration parameter:  w0 * T, (default: 10)]
%  CR_win - chirp rate of the window
%  IFE_method - instantaneous frequency estimation method
%  CRE_method - chirp rate estimation method (1-4)
%  Fourier_method - Fourier transform computation method - point by point or FFT
%  OUT:
%  S - STFT
%  HSS - horizontally synchrosqueezed STFT
%  CR - chirp rate distribution

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
if ~exist('IFE_method', 'var')
    IFE_method = 1;
end
if ~exist('CRE_method', 'var')
    CRE_method = 2;
end
if ~exist('Fourier_method', 'var')
    Fourier_method = 'FFT';
end

S = zeros(N_FFT, x.N);
SdW = zeros(N_FFT, x.N);
Sd2W = zeros(N_FFT, x.N);
StW = zeros(N_FFT, x.N);
St2W = zeros(N_FFT, x.N);
St2dW = zeros(N_FFT, x.N);
StdW = zeros(N_FFT, x.N);
Std2W = zeros(N_FFT, x.N);
CR = zeros(N_FFT, x.N);
alpha = zeros(N_FFT, x.N);
HSS = zeros(N_FFT, x.N);
if(CRE_method == 4)
    Sd3W   = zeros(N_FFT, x.N);
    St2d2W = zeros(N_FFT, x.N);
end

K = N_FFT;
f_points = 1:N_FFT;
f_points = f_points - N_FFT/2;
A = -1j * 2*pi / N_FFT;

if strcmp(Fourier_method, 'ptByPt')
    for n = 1:x.N
        k_min = min(n-1, round(K/2));
        k_max = min(x.N-n, round(K/2));
        k = ((-k_min):k_max).';
        w     = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
        dw    = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
        d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
        tw    = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
        t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
        t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR_win);
        td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR_win);
        tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
        
        if CRE_method == 4 
            d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR_win);
            t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR_win);
        end
        
        for m = 1:N_FFT
            S(m,n)      = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* w     .* exp(A .* f_points(m) .* k));
            SdW(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* dw    .* exp(A .* f_points(m) .* k));
            Sd2W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* d2w   .* exp(A .* f_points(m) .* k));
            StW(m,n)    = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tw    .* exp(A .* f_points(m) .* k));
            StdW(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* tdw   .* exp(A .* f_points(m) .* k));
            Std2W(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* td2w  .* exp(A .* f_points(m) .* k));
            St2W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2w   .* exp(A .* f_points(m) .* k));
            St2dW(m,n)  = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2dw  .* exp(A .* f_points(m) .* k));
            
            if CRE_method == 4 
                Sd3W(m,n)   = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* d3w   .* exp(A .* f_points(m) .* k));
                St2d2W(m,n) = exp(A * f_points(m) * (n-1)) * sum(signal(n+k) .* t2d2w .* exp(A .* f_points(m) .* k));
            end
            
            %% general reassignment coordinates
            t = n + StW(m,n) / S(m,n);
            f = 1j * m - N_FFT/(2*pi) * SdW(m,n) / S(m,n);
            
            %% chirp rate estimation
            if CRE_method == 1
                num = real(SdW(m,n) * conj(S(m,n)));
                den = imag(StW(m,n) * conj(S(m,n)));
                
                alpha(m,n) = 1j * num/den;
            elseif CRE_method == 2
                num = (Sd2W(m,n) * S(m,n) - SdW(m,n)^2);
                den = StW(m,n) * SdW(m,n) - StdW(m,n) * S(m,n);
                
                alpha(m,n) = num/den;
            elseif CRE_method == 3
                num = ( StdW(m,n) * S(m,n) -  StW(m,n) * SdW(m,n) + S(m,n)^2) ;
                den = StW(m,n)^2 - St2W(m,n) * S(m,n);
                
                alpha(m,n) = num/den;
            elseif CRE_method == 4
                num = -(St2W(m,n) * Sd2W(m,n)^2 - St2dW(m,n) * Sd2W(m,n) * SdW(m,n) - St2d2W(m,n) * S(m,n) * Sd2W(m,n) + St2d2W(m,n) * SdW(m,n)^2 - Sd3W(m,n) * St2W(m,n) * SdW(m,n) + Sd3W(m,n) * St2dW(m,n) * S(m,n));
                den = (Sd2W(m,n) * StdW(m,n) * St2W(m,n) - Sd2W(m,n) * St2dW(m,n) * StW(m,n) - SdW(m,n) * St2W(m,n) * Std2W(m,n) + SdW(m,n) * St2d2W(m,n) * StW(m,n) - StdW(m,n) * St2d2W(m,n) * S(m,n) + St2dW(m,n) * Std2W(m,n) * S(m,n));
                alpha(m,n) = num/den;
            else
                error('Wrong CR estimator type!')
            end
            CR(m,n) = imag(alpha(m,n));
            
             %% group delay estimation
            if IFE_method == 1 
                %% classical GD (reassignment operator)
                precise_GD  = real(round(t));
            elseif IFE_method == 2
                if CR(m,n) ~= 0
                    precise_GD = round(imag(alpha(m,n) * t) / imag(alpha(m,n)) + round((2*pi)/N_FFT * imag(alpha(m,n))^(-1) * (m - imag(f))));
                else
                    precise_GD = round(real(t));
                end
            else
                error('Wrong IF estimator type!')
            end
            if precise_GD < 1 || precise_GD > x.N || isnan(precise_GD)
                continue;
            end 
            HSS(m,precise_GD) = HSS(m,precise_GD) + S(m,n);
        end
    end
elseif strcmp(Fourier_method, 'FFT')
     k = -K/2:K/2-1;
    w     = Gab_Gaussian_Window(k, L, 0, 0, CR_win);
    dw    = Gab_Gaussian_Window(k, L, 1, 0, CR_win);
    d2w   = Gab_Gaussian_Window(k, L, 2, 0, CR_win);
    tw    = Gab_Gaussian_Window(k, L, 0, 1, CR_win);
    t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR_win);
    t2dw  = Gab_Gaussian_Window(k, L, 1, 2, CR_win);
    td2w  = Gab_Gaussian_Window(k, L, 2, 1, CR_win);
    tdw   = Gab_Gaussian_Window(k, L, 1, 1, CR_win);
    if(CRE_method == 4 || IFE_method == 3 || IFE_method == 4 || IFE_method == 5)
        d3w   = Gab_Gaussian_Window(k, L, 3, 0, CR_win);
        t2d2w = Gab_Gaussian_Window(k, L, 2, 2, CR_win);
    end
    signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
    phase_shift = linspace(1/N_FFT,K,N_FFT).';
    n_bins = 1:x.N;
    for n = 1:x.N
  tmp = signal(n_bins(n):n_bins(n)+K-1);
        S(:,n)      = fftshift(fft((tmp .* w),    N_FFT));
        SdW(:,n)    = fftshift(fft((tmp .* dw),   N_FFT));
        Sd2W(:,n)   = fftshift(fft((tmp .* d2w),  N_FFT));
        StW(:,n)    = fftshift(fft((tmp .* tw),   N_FFT));
        StdW(:,n)   = fftshift(fft((tmp .* tdw),  N_FFT));
        Std2W(:,n)  = fftshift(fft((tmp .* td2w), N_FFT));
        St2W(:,n)   = fftshift(fft((tmp .* t2w),  N_FFT));
        St2dW(:,n)  = fftshift(fft((tmp .* t2dw), N_FFT));
        
        if CRE_method == 4 
            Sd3W(:,n)   = fftshift(fft(signal(n+k) .* d3w, N_FFT));
            St2d2W(:,n) = fftshift(fft(signal(n+k) .* t2d2w, N_FFT));
        end
        
        for m = 1:N_FFT
            %% general reassignment coordinates%             
            t = n + StW(m,n) / S(m,n);
            f = 1j * m -  N_FFT/(2*pi) *SdW(m,n) / S(m,n);
                        
            %% chirp rate estimation
            if CRE_method == 1
                num = real(SdW(m,n) * conj(S(m,n)));
                den = imag(StW(m,n) * conj(S(m,n)));
                
                alpha(m,n) = 1j * num/den;
            elseif CRE_method == 2
                num = (Sd2W(m,n) * S(m,n) - SdW(m,n)^2);
                den = StW(m,n) * SdW(m,n) - StdW(m,n) * S(m,n);
                
                alpha(m,n) = num/den;
            elseif CRE_method == 3
                num = ( StdW(m,n) * S(m,n) -  StW(m,n) * SdW(m,n) + S(m,n)^2) ;
                den = StW(m,n)^2 - St2W(m,n) * S(m,n);
                
                alpha(m,n) = num/den;
            elseif CRE_method == 4
                num = -(St2W(m,n) * Sd2W(m,n)^2 - St2dW(m,n) * Sd2W(m,n) * SdW(m,n) - St2d2W(m,n) * S(m,n) * Sd2W(m,n) + St2d2W(m,n) * SdW(m,n)^2 - Sd3W(m,n) * St2W(m,n) * SdW(m,n) + Sd3W(m,n) * St2dW(m,n) * S(m,n));
                den = (Sd2W(m,n) * StdW(m,n) * St2W(m,n) - Sd2W(m,n) * St2dW(m,n) * StW(m,n) - SdW(m,n) * St2W(m,n) * Std2W(m,n) + SdW(m,n) * St2d2W(m,n) * StW(m,n) - StdW(m,n) * St2d2W(m,n) * S(m,n) + St2dW(m,n) * Std2W(m,n) * S(m,n));
                alpha(m,n) = num/den;
            else
                error('Wrong CR estimator type!')
            end
            CR(m,n) = imag(alpha(m,n));
            
            %% group delay estimation
            if IFE_method == 1 
                %% classical GD (reassignment operator)
                precise_GD  = real(round(t));
            elseif IFE_method == 2
                if CR(m,n) ~= 0
                    precise_GD = round(imag(alpha(m,n) * t) / imag(alpha(m,n)) + round((2*pi)/N_FFT * imag(alpha(m,n))^(-1) * (m - imag(f))));
                else
                    precise_GD = round(real(t));
                end
            else
                error('Wrong IF estimator type!')
            end
            if precise_GD < 1 || precise_GD > x.N || isnan(precise_GD)
                continue;
            end 
            S(m,n) = S(m,n)*exp(1j * pi * phase_shift(m));
            HSS(m,precise_GD) = HSS(m,precise_GD) + S(m,n) *exp(-1j*2*pi*f_points(m)*(n-1)/N_FFT);
        end
    end
    
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end
end



