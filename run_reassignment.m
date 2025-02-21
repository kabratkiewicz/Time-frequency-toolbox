%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., & Samczyński, P. (2021, September). 
%   Multitaper time-frequency reassigned spectrogram in micro-Doppler radar 
%   signal analysis. In 2021 Signal Processing Symposium (SPSympo) (pp. 1-5). IEEE.

close all
clear
clc

fontsize = 20;    % fontsize
img_max_size = 1; % 1 - fullscreen, otherwise - default size

addpath('UTILS/')
addpath('GAB/')
Init_Env(fontsize, img_max_size)

fs = 512;
T = 1;
t = 0:1/fs:T-1/fs;
[a, if1,if2,s1,s2,signal] = two_cos_sig(t);
signal = awgn(signal, 20, 'measured'); % Signal Processing Toolbox required

sigma = 5;
gamma_K = 1e-4;
N_FFT = 1024;
CR_win = 0;%1.5e11;%4.0e9
threshold = -50;
method = 'FFT';
fs = 1;
T = length(signal);

x.signal = signal(:);
x.fs = fs;
x.N = length(signal);
x.T = T;

%%
S = Gab_STFT(x, N_FFT, sigma, 0,method);
f_scale = linspace(-0.5, 0.5, N_FFT) ;
t_scale = 1:x.N;
Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, 1);

IFreq = Gab_Get_IFreq_Est(1, x, N_FFT , sigma, 0,method);
GDel = Gab_Get_Spectral_Delay(x, N_FFT, sigma, 0,method);
R = Gab_TF_Reassignment(GDel, IFreq, S);
Plot_Energy(R, threshold, 0, t_scale, f_scale, fontsize, 1,'power');

