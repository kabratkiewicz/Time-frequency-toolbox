function H = Plot_AS(AS, AS_min, AS_max, t_scale, f_scale, fontsize, max_size, lin_log)

if (max_size == 1)
    H = figure('units','normalized','outerposition',[0 0 1 1]);
else
    H = figure;
end

imagesc(t_scale,f_scale,AS./1e3);
set(gca,'YDir','normal')
xlabel('t [s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [Hz]','FontSize',fontsize, 'interpreter','latex')
h=colorbar;title(h,'$\hat{\gamma}_x(t,\omega)$ [kHz/$^3$s]','FontSize',fontsize, 'interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')
clim([AS_min AS_max])
colormap(jet)
c = colorbar;
c.Label.String = '$\hat{\gamma}_x(t,\omega)$ [kHz/s$^3$]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
if strcmp(lin_log,'log')
    set(gca,'ColorScale','log')
else
    set(gca,'ColorScale','lin')
end

drawnow
end