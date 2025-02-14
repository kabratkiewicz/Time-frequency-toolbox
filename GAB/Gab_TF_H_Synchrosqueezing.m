function SE = Gab_TF_H_Synchrosqueezing(S, sDelay, method)
%  S - STFT
%  sDelay - instantaneous frequency estimate
%  method - Fourier transform computation method - point by point or FFT
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Abratkiewicz, K., Gambrych, J., Stasiak, K., & Samczyński, P. 
%   (2021, September). Estimation of Rotational Speed of Helicopter Rotor Through 
%   Horizontal Synchrosqueezing. In 2021 Signal Processing Symposium (SPSympo) (pp. 6-10). IEEE.

if size(S) ~= size(sDelay)
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
            t_idx = sDelay(j,i);
            if t_idx < 1 || t_idx > N
                continue;
            end
            SE(j,t_idx) = SE(j,t_idx) + S(j,i);%/(2*pi) * exp(1j*2*pi*m(j)*(i-1)/M);
        end
    end
elseif strcmp(method, 'FFT')
    for i = 1:N
        for j = 1:M
            t_idx = sDelay(j,i);
            if t_idx < 1 || t_idx > N
                continue;
            end
            SE(j,t_idx) = SE(j,t_idx) + S(j,i);
        end
    end
else
    error('Wrong transform method. Available methods ptByPt or FFT');
end

end