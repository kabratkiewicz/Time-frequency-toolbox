function H = Plot_Energy(S_energy, e_min, e_max, t_scale, f_scale, fontsize, max_size, sig_type)

if ~exist('sig_type', 'var')
    sig_type = 'voltage';
end

S = db(abs(S_energy), sig_type) - max(max(db(abs(S_energy), sig_type)));

if (max_size == 1)
    H = figure('units','normalized','outerposition',[0 0 1 1]);
else
    H = figure;
end

imagesc(t_scale, f_scale, S,'interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('Time samples','FontSize',fontsize, 'interpreter','latex')
ylabel('Normalized freq.','FontSize',fontsize, 'interpreter','latex')
h=colorbar;title(h,'E[dB]','FontSize',fontsize, 'interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')
colormap(flipud(gray))
clim([e_min e_max]) 
c = colorbar;
c.Label.String = 'E[dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
drawnow
end

