%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., & Gambrych, J. (2022). 
%   Real-time variants of vertical synchrosqueezing: Application 
%   to radar remote sensing. 
%   IEEE Journal of Selected Topics in Applied Earth Observations 
%   and Remote Sensing, 15, 1760-1774.

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
signal = awgn(signal, 50, 'measured'); % Signal Processing Toolbox required

sigma = 5;
N_FFT = 1024;
CR_win = 0;
threshold = -50;
method = 'FFT';
T = length(signal);

clear x;

x.signal = signal(:);
x.fs = fs;
x.N = length(signal);
x.T = T;

%%
S = Gab_STFT(x, N_FFT, sigma, 0, method);
f_scale = linspace(-0.5, 0.5, N_FFT) ;
t_scale = 1:x.N;
Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(S).^2)
x_rec_S = Gab_ISSTFT(S, sigma, N_FFT);
%% 1st-order vertical synchrosqueezing
IFreq1 = Gab_Get_IFreq_Est(1, x, N_FFT , sigma, 0,method);
VSS1 = Gab_TF_V_Synchrosqueezing(S,IFreq1,method);
Plot_Energy(VSS1, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(VSS1).^2)
x_rec_VSS1 = Gab_ISSTFT(S, sigma, N_FFT);
%% 2nd-order vertical synchrosqueezing
IFreq2 = Gab_Get_IFreq_Est(2, x, N_FFT , sigma, 0,method);
VSS2 = Gab_TF_V_Synchrosqueezing(S,IFreq2,method);
Plot_Energy(VSS2, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(VSS2).^2)
x_rec_VSS2 = Gab_ISSTFT(S, sigma, N_FFT);
%% 3rd-order vertical synchrosqueezing
IFreq3 = Gab_Get_IFreq_Est(3, x, N_FFT , sigma, 0,method);
VSS3 = Gab_TF_V_Synchrosqueezing(S,IFreq3,method);
Plot_Energy(VSS3, threshold, 0, t_scale, f_scale, fontsize, 1);
Renyi_Entropy(abs(VSS3).^2)
x_rec_VSS3 = Gab_ISSTFT(S, sigma, N_FFT);