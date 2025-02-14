%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., & Samczyński, P. (2021, September). 
%   Multitaper time-frequency reassigned spectrogram in micro-Doppler radar 
%   signal analysis. In 2021 Signal Processing Symposium (SPSympo) (pp. 1-5). IEEE.

close all
clear
clc

addpath("GAB\")
addpath("UTILS\")
addpath("MULTAP\")

threshold = -50;
fontsize = 20;
Init_Env(fontsize,1)

fs = 512;
T = 1;
t = 0:1/fs:T-1/fs;
[a, if1,if2,s1,s2,signal] = two_cos_sig_slow(t);
signal = awgn(signal, 20, 'measured'); % Signal Processing Toolbox required

x.signal = signal;
x.N = length(signal);
x.fs = 1;
N_FFT = 1024;
sigma = 80;
M = 4;

f_scale = linspace(-x.fs/2, x.fs/2, N_FFT) ;
t_scale = linspace(0,x.N/x.fs,x.N);

[S, iFreq, gDel] = Gab_Multitaper_Reass_Vectors(x, N_FFT, sigma, M);

R = zeros(N_FFT, x.N, M);

for i = 1:M
    R(:,:,i) = Gab_TF_Reassignment(gDel(:,:,i), iFreq(:,:,i), S(:,:,i));
    Plot_Energy(R(:,:,i), threshold, 0, t_scale, f_scale, fontsize, 1,'power');
end

R_mean = mean(R,3);
S_mean = mean(abs(S).^2,3);
Plot_Energy(S_mean, threshold, 0, t_scale, f_scale, fontsize, 1,'power');
Plot_Energy(R_mean, threshold, 0, t_scale, f_scale, fontsize, 1,'power');