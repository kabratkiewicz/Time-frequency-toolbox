%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: K. Abratkiewicz, B. Falęcki, A. Burzyńska, M. Pożoga, H. Rothkaehl
%   "On the Influence of Ionospheric-Induced Polarization Distortion on
%   Over-the-Horizon Waves" International Radar Symposium 2026

close all
clear
clc

addpath("GAB\")
addpath("UTILS\")
load("colormap_BR.mat")

NFFT = 512;

gamma = 1e-4;
method = 'FFT';
L = 10;
img_max_size = 1;
fontsize = 30;
threshold = -50;
Init_Env(fontsize, img_max_size);

N = 1024;
t = linspace(0,1,N);

phi1 = -700 * pi .* t + 700 * pi .* t.^2;
a1 = 1;
xH = a1 .* exp(1i * phi1);
xH = awgn(xH, 50, "measured");

phi2 = -700 * pi .* t + 700 * pi .* t.^2 + (5 .*t * pi);
a2 = 1;
xV = a2 .* exp(1i * phi2);
xV = awgn(xV, 50, "measured");

x.signal = xH;
x.N = length(xH);
x.fs = N;
SH = Gab_STFT(x, NFFT, L, 0, method);
f_scale = linspace(-x.fs/2, x.fs/2, NFFT);

f_scale = f_scale./1e3;
t_scale = linspace(0,x.N/x.fs,x.N);

x.signal = xV;
x.N = length(xV);
x.fs = N;
SV = Gab_STFT(x, NFFT, L, 0, method);

[S0, S1n, S2n, S3n] = Stokes_Parameter_estimation(SH, SV);

figure;
imagesc(t_scale, f_scale, S0)
set(gca,'ydir','normal')
xlabel('t [s]','FontSize', fontsize, 'interpreter','latex')
ylabel('f [kHz]','FontSize', fontsize, 'interpreter','latex')

colormap('turbo')
col=colorbar;
col.Label.Interpreter = 'latex';
col.TickLabelInterpreter = 'latex';
col.Label.String = '$S_0$';
max_energy = (max(max(S0)));
d = db(S0/max_energy);
q = (threshold - d) / (threshold);
clim([0 1])

figure;
im=imagesc(t_scale, f_scale, S1n);
set(gca,'ydir','normal')
xlabel('t [s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [kHz]','FontSize',fontsize, 'interpreter','latex')
colormap(cmap)
col=colorbar;
col.Label.Interpreter = 'latex';
col.TickLabelInterpreter = 'latex';
col.Label.String = '$S_1$';
im.AlphaData = q;
clim([-1 1])

figure;
im=imagesc(t_scale, f_scale, S2n);
set(gca,'ydir','normal')
xlabel('t [s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [kHz]','FontSize',fontsize, 'interpreter','latex')
colormap(cmap)
col=colorbar;
col.Label.Interpreter = 'latex';
col.TickLabelInterpreter = 'latex';
col.Label.String = '$S_2$';
im.AlphaData = q;
clim([-1 1])

figure;
im=imagesc(t_scale, f_scale,S3n);
set(gca,'ydir','normal')
xlabel('t [s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [kHz]','FontSize',fontsize, 'interpreter','latex')
colormap(cmap)
col=colorbar;
col.Label.Interpreter = 'latex';
col.TickLabelInterpreter = 'latex';
col.Label.String = '$S_3$';
im.AlphaData = q;
clim([-1 1])

