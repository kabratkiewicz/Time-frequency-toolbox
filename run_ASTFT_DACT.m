%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K. (2020), Double-adaptive chirplet transform for 
%   radar signature extraction. IET Radar Sonar Navig., 14: 1463-1474. 
%   https://doi.org/10.1049/iet-rsn.2020.0084

close all
clear
clc

addpath("GAB\")
addpath("UTILS\")

threshold = -50;
fontsize = 50;
Init_Env(50,1)

fs = 1000;
N = 1;
t = 0:1/fs:1-1/fs;
t = t - mean(t);
signal = exp(1j*2*pi*250.*t.^2+ 1j*2*pi*200.*t);
signal = awgn(signal,20,'measured');

x.signal = signal;
x.N = length(signal);
x.fs = fs;
method = 'FFT';
NFFT = 1024;
L = 10;
S = Gab_STFT(x, NFFT, L, 0, method);
CR_est = 3;

f_scale = linspace(-0.5, 0.5, N_FFT) ;
t_scale = 1:x.N;
ff = Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(S).^2)
%% chirp rate estimation
CR = Gab_Calculate_CR(x, NFFT, L, 0, 0, CR_est, method);
[~,ICR] = Get_ICR_Vector(S,CR,10);
ASTFT_sigma = sqrt(1./2./pi./(abs(ICR)));
%% ASTFT
AS = Gab_ASTFT(x, NFFT, ASTFT_sigma, 0, method);
Plot_Energy(AS, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(AS).^2)
%% Chirplet transform
CT = Gab_CT(x, NFFT, L, 0, ICR, method);
Plot_Energy(CT, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(CT).^2)
%% Double adaptive chirplet transform
DACT = Gab_ACT(x, NFFT, ASTFT_sigma, 0, ICR, method);
Plot_Energy(DACT, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(DACT).^2)
