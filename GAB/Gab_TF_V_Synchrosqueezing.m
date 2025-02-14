function SE = Gab_TF_V_Synchrosqueezing(S, IFreq, method)
%  S - STFT
%  IFreq - instantaneous frequency estimate
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Drozdowicz, J., & Samczyński, P. (2023). 
%   Vertical Synchrosqueezing for High-Resolution Radar Imaging. 
%   IEEE Transactions on Geoscience and Remote Sensing.

if size(S) ~= size(IFreq)
    error('Matrices have to be the same size')
end
if ~exist('method', 'var')
    method = 'FFT';
end
SE = zeros(size(S));
N = size(S,2); % time samples
M = size(S,1); % frequency bins
m = 1:M;
m = m - M/2;
if strcmp(method, 'ptByPt')
    for i = 1:N
        for j = 1:M
            f_idx = IFreq(j,i);
            if f_idx < 1 || f_idx > M
                continue;
            end
            SE(f_idx,i) = SE(f_idx,i) + S(j,i) * exp(1j*2*pi*m(j)*(i-1)/M);
        end
    end
elseif strcmp(method, 'FFT')
    for i = 1:N
        for j = 1:M
            f_idx = IFreq(j,i) + 1;
            if f_idx < 1 || f_idx > M
                continue;
            end
            SE(f_idx,i) = SE(f_idx,i) + S(j,i);
        end
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end
end

