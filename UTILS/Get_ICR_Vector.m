function [idx,ICR] = Get_ICR_Vector(S,CR, s_order)
% S - spectrogram
% CR - accelerogram
    [~,idx] = max(S,[],1);
    ICR = zeros(1,length(idx));
    for i=1:length(idx)
        ICR(i) = CR(idx(i), i);
    end
    
    for i=1:s_order
        ICR = smooth(medfilt1(ICR));
    end
end

