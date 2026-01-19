function [iFreq, alpha, beta, gamma] = Gab_Get_4th_Order_Param_Est(x, N_FFT, L, CR)
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
%   related paper: Abratkiewicz, K., Drozdowicz, J., & Samczyński, P. (2023). 
%   Vertical Synchrosqueezing for High-Resolution Radar Imaging. 
%   IEEE Transactions on Geoscience and Remote Sensing.

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

S        = zeros(N_FFT, x.N);
StW      = zeros(N_FFT, x.N);
SdW      = zeros(N_FFT, x.N);
St2W     = zeros(N_FFT, x.N);
St3W     = zeros(N_FFT, x.N);
St4W     = zeros(N_FFT, x.N);
St5W     = zeros(N_FFT, x.N);
St6W     = zeros(N_FFT, x.N);
SdtW     = zeros(N_FFT, x.N);
Sdt2W     = zeros(N_FFT, x.N);
Sdt3W     = zeros(N_FFT, x.N);

iFreq = zeros(N_FFT, x.N);
alpha = zeros(N_FFT, x.N);
beta  = zeros(N_FFT, x.N);
gamma = zeros(N_FFT, x.N);

K = N_FFT;
k = -K/2:K/2-1;

w     = Gab_Gaussian_Window(k, L, 0, 0, CR);
tw    = Gab_Gaussian_Window(k, L, 0, 1, CR);
t2w   = Gab_Gaussian_Window(k, L, 0, 2, CR);
t3w   = Gab_Gaussian_Window(k, L, 0, 3, CR);
t4w   = Gab_Gaussian_Window(k, L, 0, 4, CR);
t5w   = Gab_Gaussian_Window(k, L, 0, 5, CR);
t6w   = Gab_Gaussian_Window(k, L, 0, 6, CR);
dtw   = Gab_Gaussian_Window(k, L, 1, 1, CR);
dt2w   = Gab_Gaussian_Window(k, L, 1, 2, CR);
dt3w   = Gab_Gaussian_Window(k, L, 1, 3, CR);
dw    = Gab_Gaussian_Window(k, L, 1, 0, CR);

signal  =  [zeros(K/2, 1); signal; zeros(K/2, 1)];
phase_shift = linspace(0,K-1,N_FFT).';
n_bins = 1:x.N;

for n = 1:x.N
    tmp = signal(n_bins(n):n_bins(n)+K-1);
    S(:,n)       = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* w,     N_FFT));
    StW(:,n)     = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* tw,    N_FFT));
    St2W(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t2w,   N_FFT));
    St3W(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t3w,   N_FFT));
    St4W(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t4w,   N_FFT));
    St5W(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t5w,   N_FFT));
    St6W(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* t6w,   N_FFT));
    SdW(:,n)     = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dw,    N_FFT));
    SdtW(:,n)    = exp(1i * pi * phase_shift) .* fftshift(fft(tmp .* dtw,   N_FFT));       
    Sdt2W(:,n)    = exp(1i * pi * phase_shift) .*fftshift( fft(tmp .* dt2w,   N_FFT));     
    Sdt3W(:,n)    = exp(1i * pi * phase_shift) .*fftshift( fft(tmp .* dt3w,   N_FFT));     
    

    for m = 1:N_FFT
        if abs(S(m,n)) > eps
            Am = [St3W(m,n)   St2W(m,n)   StW(m,n)   S(m,n);
                St4W(m,n)   St3W(m,n)   St2W(m,n)   StW(m,n);
                St5W(m,n)   St4W(m,n)   St3W(m,n)   St2W(m,n);
                St6W(m,n)   St5W(m,n)   St4W(m,n)   St3W(m,n);];

            Bm = [ SdW(m,n);
                SdtW(m,n) + S(m,n);
                Sdt2W(m,n) + 2*StW(m,n);
                Sdt3W(m,n) + 3*St2W(m,n)];

            tau = (StW(m,n)/S(m,n));

            Cm = Am \ Bm;
            hat_v3 = -Cm(1);
            hat_v2 = -3*hat_v3 * tau  - Cm(2);
            hat_v1 = -3*hat_v3*tau^2 - 2*hat_v2*tau - Cm(3);
            hat_v0 = -hat_v3 * tau^3 - hat_v2*tau^2 - hat_v1*tau - Cm(4);
            gamma(m,n) = imag(hat_v3)/2/pi;
            beta(m,n) = imag(hat_v2)/2/pi;
            alpha(m,n) = imag(hat_v1)/2/pi;
            iFreq(m,n) = m+round(N_FFT*imag(hat_v0)/2/pi);

        end
    end
end

end


