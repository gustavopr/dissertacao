%% GRÁFICOS DA FUNCAO PERDA ASSIMETRICA

lambda = 1;
gamma = -1.0;
alfa = 1.0;

h = figure;

% gráfico da perda em relação a desvios do produto

l_hiato = @(x) lambda * (exp(gamma * x) - gamma * x - 1)./gamma^2;

subplot(1,2,1);
fplot(l_hiato,[-1.5,2.5],'color','black');
grid
axis tight
xlabel('(a)'); % x-axis label
ylabel('Perda'); % y-axis label
legend('gamma = -1.0');

% gráfico da perda em relação a desvios da inflação à meta

l_hiato = @(pi) (exp(alfa * pi) - alfa * pi - 1)./alfa^2;

subplot(1,2,2);
fplot(l_hiato,[-2.5,1.5],'color','black');
grid
axis tight
xlabel('(b)'); % x-axis label
ylabel('Perda'); % y-axis label
legend('alfa = 1.0');

set(h, 'Position', [0 0 900 300])
h.PaperPosition = [0 0 900 300];
h.PaperPositionMode = 'manual';

%% GRÁFICOS DA FUNCAO LOGISTICA PARA DIFERENTES GAMMAS



h = figure(2);
% gráfico da funcao exponencial com Gamma positivo e suave

c = 0;
Gamma = 1;

l_exp = @(z) 1/(1+exp(Gamma*(z - c)));

subplot(2,2,1);
fplot(l_exp,[-10,10],'color','black');
ylim([-0.2 1.2]);
grid
xlabel('(a)'); % x-axis label
legend('Gamma = 1.0', 'Location', 'NorthEast');

% gráfico da funcao exponencial com Gamma negativo e suave

Gamma = -1;

l_exp = @(z) 1/(1+exp(Gamma*(z - c)));

subplot(2,2,2);
fplot(l_exp,[-10,10],'color','black');
ylim([-0.2 1.2]);
grid
xlabel('(b)'); % x-axis label
legend('Gamma = -1.0', 'Location', 'SouthEast');

% gráfico da funcao exponencial com Gamma positivo e  não suave

Gamma = 10;

l_exp = @(z) 1/(1+exp(Gamma*(z - c)));

subplot(2,2,3);
fplot(l_exp,[-10,10],'color','black');
ylim([-0.2 1.2]);
grid
xlabel('(c)'); % x-axis label
legend('Gamma = 10.0', 'Location', 'NorthEast');

% gráfico da funcao exponencial com Gamma negativo e  não suave

Gamma = -10;

l_exp = @(z) 1/(1+exp(Gamma*(z - c)));

subplot(2,2,4);
fplot(l_exp,[-10,10],'color','black');
ylim([-0.2 1.2]);
grid
xlabel('(d)'); % x-axis label
legend('Gamma = -10.0', 'Location', 'SouthEast');

set(h, 'Position', [0 0 900 700])
h.PaperPosition = [0 0 900 700];
h.PaperPositionMode = 'manual';


% salvando
%saveas(h, 'D:\Documents\Dissertacao\imagens\funcao_perda_assimetrica.png','png');
%print('D:\Documents\Dissertacao\imagens\func_exponencial_gamma_pos_neg.png','-dpng','-r0')

%% GRÁFICOS DAS VARIAVEIS UTILIZADAS NO MODELO

data = x2mdate(x(:,1));

h = figure(3);

subplot(2,2,1);
plot(data,x(:,4),'-k');
xlabel('Inflação esperada'); % x-axis label
ylabel('%'); % y-axis label
grid
xlim([731582 736024]);
dateaxis('x',12)

subplot(2,2,2);
plot(data,x(:,7),'-k');
xlabel('Diferencial de inflação à meta'); % x-axis label
ylabel('%'); % y-axis label
grid
xlim([731582 736024]);
dateaxis('x',12)

subplot(2,2,3);
plot(data,x(:,3),'-k');
xlabel('Taxa Selic Efetiva'); % x-axis label
ylabel('%'); % y-axis label
grid
xlim([731582 736024]);
dateaxis('x',12)

subplot(2,2,4);
plot(data,x(:,2),'-k');
xlabel('Hiato do Produto'); % x-axis label
grid
xlim([731582 736024]);
dateaxis('x',12)

set(h, 'Position', [0 0 900 700])
h.PaperPosition = [0 0 900 700];
h.PaperPositionMode = 'manual';
