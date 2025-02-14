function [profile] = Calculate_R_profile(s_energy, c_rate, start, stop, points)
    
  scale = linspace(start,stop,points+1);
 
  profile = zeros(1,points);
  for i = 1:size(s_energy,2) % kolumna - t
        for j = 1:size(s_energy,1) % wiersz - f
            cr = c_rate(j,i);
            for k = 1:length(scale)-1
                if (cr > scale(k)) && (cr < scale(k+1))
                    profile(k) = profile(k) + s_energy(j,i);
                end
            end
        end
        
  end

  
end

