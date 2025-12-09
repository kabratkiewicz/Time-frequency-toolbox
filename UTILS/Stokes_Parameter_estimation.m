function [S0, S1n, S2n, S3n] = Stokes_Parameter_estimation(SH,SV)
% STOKES_PARAMETER_ESTIMATION  Compute Stokes parameters from dual-polarization signals
%
%   [S0, S1n, S2n, S3n] = STOKES_PARAMETER_ESTIMATION(SH, SV)
%   estimates the four Stokes parameters based on two complex input signals
%   representing horizontal (SH) and vertical (SV) polarizations.
%
%   The function computes:
%       S0 – total intensity of the electromagnetic wave,
%       S1 – difference between horizontal and vertical power,
%       S2 – correlation term describing linear polarization at ±45°,
%       S3 – correlation term describing circular polarization.
%
%   The parameters S1, S2, and S3 are then normalized by S0 (yielding S1n,
%   S2n, and S3n). Additionally, S0 is scaled by its maximum value for
%   convenient visualization.
%
%   INPUTS:
%       SH  - complex signal in horizontal polarization
%       SV  - complex signal in vertical polarization
%
%   OUTPUTS:
%       S0  - normalized total power
%       S1n - normalized S1 parameter
%       S2n - normalized S2 parameter
%       S3n - normalized S3 parameter
%
%   AUTHOR:
%       Karol Abratkiewicz
%       email: karol.abratkiewicz@pw.edu.pl
%
%   REFERENCE:
%       K. Abratkiewicz, B. Falęcki, A. Burzyńska, M. Pożoga, H. Rothkaehl
%       "On the Influence of Ionospheric-Induced Polarization Distortion on
%       Over-the-Horizon Waves" IEEE Antennas and Wireless Propagation Letters
%       
%
%   ---------------------------------------------------------------------

S0 = abs(SH).^2 + abs(SV).^2;
S1 = abs(SH).^2 - abs(SV).^2;
S2 = 2 * real(SH .* conj(SV));
S3 = 2 * imag(SH .* conj(SV));

S1n = S1 ./ S0;
S2n = S2 ./ S0;
S3n = S3 ./ S0;
S0 = S0./max(max(S0));
end