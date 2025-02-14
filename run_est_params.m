%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Drozdowicz, J., & Samczyński, P. (2023). 
%   Vertical Synchrosqueezing for High-Resolution Radar Imaging. 
%   IEEE Transactions on Geoscience and Remote Sensing.

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
[a, if1,if2,s1,s2,signal] = two_cos_sig_slow(t);
signal = awgn(signal, 100);

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

%%
S = Gab_STFT(x, N_FFT, sigma, 0,method);
f_scale = linspace(-x.fs/2, x.fs/2, N_FFT);
t_scale = linspace(0,x.N/x.fs,x.N);
Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, 1);

%%
[~, alpha, beta, gamma] = Gab_Get_4th_Order_Param_Est(x, N_FFT, sigma, 0);
%%
Plot_CR(alpha*fs^2, -1, 1, t_scale, f_scale, fontsize, img_max_size, 'lin')
Plot_AJ(beta*fs^3 ,-5, 5, t_scale, f_scale, fontsize, img_max_size, 'lin')
Plot_AS(gamma*fs^4,-20, 20, t_scale, f_scale, fontsize, img_max_size, 'lin')
