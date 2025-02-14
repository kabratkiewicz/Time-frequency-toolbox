function H = Plot_Saturated_CR(CR, E, CR_min, CR_max, threshold,t_scale, f_scale, fontsize, max_size, linlog)

if (max_size == 1)
    H = figure('units','normalized','outerposition',[0 0 1 1]);
else
    H = figure;
end

E = (abs(E));
max_energy = (max(max(E)));

d = db(E/max_energy);
q = (threshold - d) / (threshold);
f_scale = f_scale./1e6;
t_scale = t_scale./1e-6;

% f_scale = linspace(-1, 1, size(E,1));
% t_scale = linspace(0,length(t_scale)-1,length(t_scale));
CR = CR./1e9;
im = imagesc(t_scale,f_scale,CR);
set(gca,'YDir','normal')
xlabel('t[$\mu$s]','FontSize',fontsize, 'interpreter','latex')
ylabel('f[MHz]','FontSize',fontsize, 'interpreter','latex')
% title('Akcelerogram sygnalu','fontsize',fontsize, 'interpreter','latex');
h=colorbar;title(h,'CR[Hz/s]','FontSize',fontsize, 'interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')
caxis([CR_min./1e9 CR_max./1e9])
colormap jet
c = colorbar;
c.Label.String = 'CR[GHz/s]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
%c.Ruler.Exponent = 3;
if strcmp(linlog,'log')
    set(gca,'ColorScale', 'log');
end
im.AlphaData = q;
drawnow
end


