%%  Author: Karol Abratkewicz & Bartosz Falęcki
%   e-mail: karol.abratkiewicz@pw.edu.pl, bartosz.falecki@pw.edu.pl
%   related paper: K. Abratkiewicz, B. Falęcki, A. Burzyńska, M. Pożoga, H. Rothkaehl
%   "On the Influence of Ionospheric-Induced Polarization Distortion on
%   Over-the-Horizon Waves" IEEE Antennas and Wireless Propagation Letters

close all
clear
clc

addpath("GAB\")
addpath("UTILS\")

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

% Synthetic circular component 1
a1 = 1; % Amplitude
phi1H = -700 * pi .* t + 700 * pi .* t.^2; % Phase of horizontal polarization
phi1V = -700 * pi .* t + 700 * pi .* t.^2 + pi/2; % Phase of vertical polarization
x1H = a1.*exp(1i * phi1H); % Horizontal polarization signal (H)
x1H = awgn(x1H, 50, "measured"); % Added Gaussian noise
x1V = a1.*exp(1i * phi1V); % Vertical polarization signal (V)
x1V = awgn(x1V, 50, "measured");
x1 = x1H + x1V;

% Synthetic circular component 2 (opposite handedness)
a2 = 0.5; % Artificial ellipticity introduced
deltaPhi = 5*pi*t; % Phase velocity drift - artificial polarization rotation along chirp
phi2H = -700 * pi .* t + 700 * pi .* t.^2 + pi/2 + deltaPhi; % Circular rotation in opposite direction
phi2V = -700 * pi .* t + 700 * pi .* t.^2 + deltaPhi;
x2H = a2.*exp(1i * phi2H);
x2H = awgn(x2H, 50, "measured");
x2V = a2.*exp(1i * phi2V);
x2V = awgn(x2V, 50, "measured");
x2 = x2H + x2V;

% Separately received H and V channels
xH = x1H + x2H; % x1H and x2H superposition
xV = x1V + x2V; % x1V and x2V superposition

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
colormap('turbo')
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
colormap('turbo')
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
colormap('turbo')
col=colorbar;
col.Label.Interpreter = 'latex';
col.TickLabelInterpreter = 'latex';
col.Label.String = '$S_3$';
im.AlphaData = q;
clim([-1 1])





