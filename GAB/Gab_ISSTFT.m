function [ s_hat ] = Gab_ISSTFT(S, L, N_FFT)
%  S - STFT
%  L - window duration parameter
%  N_FFT - number of frequency bins
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Drozdowicz, J., & Samczyński, P. (2023). 
%   Vertical Synchrosqueezing for High-Resolution Radar Imaging. 
%   IEEE Transactions on Geoscience and Remote Sensing.

s_hat = sum(S,1) * sqrt(2*pi) * L/N_FFT;

end