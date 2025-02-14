%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Gambrych, J., Stasiak, K., & Samczyński, P. 
%   (2021, September). Estimation of Rotational Speed of Helicopter Rotor Through 
%   Horizontal Synchrosqueezing. In 2021 Signal Processing Symposium (SPSympo) (pp. 6-10). IEEE.

close all
clear 
clc

fontsize = 20;      
img_max_size = 1;   
threshold = -50;    

addpath('GAB/')
addpath('UTILS/')
Init_Env(fontsize, img_max_size); 

fs = 512;
T = 1;
t = 0:1/fs:T-1/fs;
[a, if1,if2,s1,s2,signal] = two_cos_sig(t);
signal = awgn(signal, 20, 'measured'); % Signal Processing Toolbox required

x.signal = signal;      
x.fs = fs;              
x.N = length(signal);   
x.T = T;            

sigma = 5;         
gamma_K = 1e-4;     
N_FFT = 1024;       
CR_win = 0;         
method = 'FFT';     

f_scale = linspace(-0.5, 0.5, N_FFT);   
t_scale = 1:x.N;     
S = Gab_STFT(x, N_FFT, sigma, CR_win*(x.fs^2), method);
Plot_Energy(S, threshold, 0, t_scale, f_scale, fontsize, img_max_size);         

GDel = Gab_Get_Spectral_Delay(x, N_FFT, sigma, CR_win*(x.fs^2), method);
HSS = Gab_TF_H_Synchrosqueezing(S, GDel, method);
Plot_Energy(HSS, threshold, 0, t_scale, f_scale, fontsize, img_max_size);  
Renyi_Entropy(abs(HSS).^2)

[~, OVS, CR] = Gab_TF_H2_Synchrosqueezing(x, N_FFT, sigma, CR_win*(x.fs^2), 2, 2, method);

Plot_Energy(OVS, threshold, 0, t_scale, f_scale, fontsize, img_max_size);
Renyi_Entropy(abs(OVS).^2)




