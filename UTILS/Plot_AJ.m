function H = Plot_AJ(AJ, AJ_min, AJ_max, t_scale, f_scale, fontsize, max_size, lin_log)

if (max_size == 1)
    H = figure('units','normalized','outerposition',[0 0 1 1]);
else
    H = figure;
end

imagesc(t_scale,f_scale,AJ./1e3);
set(gca,'YDir','normal')
xlabel('t [s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f [Hz]','FontSize',fontsize, 'interpreter','latex')
h=colorbar;title(h,'$\hat{\beta}_x(t,\omega)$ [kHz/$^2$s]','FontSize',fontsize, 'interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')
clim([AJ_min AJ_max])
colormap(jet)
c = colorbar;
c.Label.String = '$\hat{\beta}_x(t,\omega)$ [kHz/s$^2$]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
if strcmp(lin_log,'log')
    set(gca,'ColorScale','log')
else
    set(gca,'ColorScale','lin')
end
drawnow
end