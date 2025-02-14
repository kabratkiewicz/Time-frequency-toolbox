function IFreq = Gab_Get_IFreq_Est(est, x, N_FFT, L, CR, method)
    if est == 1
        IFreq = Gab_Get_IFreq_by_1st_Est(x, N_FFT, L, CR, method);
    elseif est == 2
        IFreq = Gab_Get_IFreq_by_2nd_Est(x, N_FFT, L, CR, method);
    elseif est == 3
        IFreq = Gab_Get_IFreq_by_3rd_Est(x, N_FFT, L, CR, method);
    end
end

