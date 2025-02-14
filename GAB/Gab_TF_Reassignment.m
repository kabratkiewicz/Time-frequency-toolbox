function RE = Gab_TF_Reassignment(sDelay, IFreq, S)
%  S - STFT
%  sDelay - instantaneous frequency estimate
%  IFreq - instantaneous frequency estimate
%%  Author: Karol Abratkewicz
%   e-mail: karol.abratkiewicz@pw.edu.pl
%   related paper: Płotka, M., Abratkiewicz, K., Malanowski, M., Samczyński, P., & Kulpa, K. 
%   (2020). The use of the reassignment technique in the time-frequency analysis 
%   applied in VHF-based passive forward scattering radar. Sensors, 20(12), 3434.

if size(S) ~= size(IFreq)
    error('Matrices have to be the same size')
end
RE = zeros(size(S));
N = size(S,2); % time samples
M = size(S,1); % frequency bins
for i = 1:N
    for j = 1:M
        f_idx = IFreq(j,i);
        t_idx = sDelay(j,i);
        if f_idx < 1 || f_idx > M
            continue;
        end
        if t_idx < 1 || t_idx > N
            continue;
        end
        RE(f_idx,t_idx) = RE(f_idx,t_idx) + abs(S(j,i))^2;
    end
end
end

