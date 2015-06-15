lambda = 1;
gamma = -1.0;
alfa = 1.0;

% l = lambda * (exp(gamma * x) - gamma * x - 1)./gamma^2 + (exp(alfa * pi) - alfa * pi - 1)./alfa^2;

h = figure;

% gráfico da perda em relação a desvios do produto

l_hiato = @(x) lambda * (exp(gamma * x) - gamma * x - 1)./gamma^2;

subplot(1,2,1);
fplot(l_hiato,[-1.5,2.5],'color','black');
%n=get(gca,'xtick');
%set(gca,'xticklabel',sprintf('%.1f |',n'));
%n=get(gca,'ytick');
%set(gca,'yticklabel',sprintf('%.1f |',n'));

xlabel('(a)'); % x-axis label
ylabel('Perda'); % y-axis label

legend('gamma = -1.0');

% gráfico da perda em relação a desvios da inflação à meta

l_hiato = @(pi) (exp(alfa * pi) - alfa * pi - 1)./alfa^2;

subplot(1,2,2);
fplot(l_hiato,[-2.5,1.5],'color','black');
%n=get(gca,'xtick');
%set(gca,'xticklabel',sprintf('%.1f |',n'));
%n=get(gca,'ytick');
%set(gca,'yticklabel',sprintf('%.1f |',n'));

xlabel('(b)'); % x-axis label
ylabel('Perda'); % y-axis label

legend('alfa = 1.0');

set(h, 'Position', [0 0 900 300])

h.PaperPosition = [0 0 900 300];
h.PaperPositionMode = 'manual';

% salvando
%saveas(h, 'D:\Documents\Dissertacao\imagens\funcao_perda_assimetrica.png','png');
print('D:\Documents\Dissertacao\imagens\funcao_perda_assimetrica.png','-dpng','-r0')