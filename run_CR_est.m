%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Samczyński, P., & Czarnecki, K. (2019). 
%   Radar signal parameters estimation using phase accelerogram in the 
%   time-frequency domain. IEEE Sensors Journal, 19(13), 5078-5085.

close all
clear
clc

fontsize = 20;    % fontsize
img_max_size = 1; % 1 - fullscreen, otherwise - default size

addpath('UTILS/')
addpath('GAB/')
addpath("MEXTR\")
Init_Env(fontsize, img_max_size)


fs = 512;
T = 1;
t = 0:1/fs:T-1/fs;
[a, if1,if2,s1,s2,signal] = two_chirps_sig(t);
signal = awgn(signal, 20);

sigma = 8;
N_FFT = 1024;
CR_win = 0;
threshold = -50;
method = 'FFT';
T = length(signal);

x.signal = signal(:);
x.fs = fs;
x.N = length(signal);
x.T = T;

est = 3;

start_point = 0;
stop_point = 1000;
N_points = 100;
cr_scale = linspace(start_point, stop_point, N_points);

%%
S = Gab_STFT(x, N_FFT, sigma, 0,method);
f_scale = linspace(-x.fs/2, x.fs/2, N_FFT);
t_scale = linspace(0,x.N/x.fs,x.N);
Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, 1);

%%
alpha = Gab_Calculate_CR(x, N_FFT, sigma, CR_win, est, method);
%%
Plot_CR(alpha*fs^2, -1, 1, t_scale, f_scale, fontsize, img_max_size, 'lin');
profile = Calculate_R_profile(abs(S), alpha*fs^2, start_point, stop_point, N_points);

figure;
plot(cr_scale, profile, 'k', 'LineWidth',2)
grid on
grid minor
xlabel('Chirp rate [Hz/s]')
ylabel('Energy')
