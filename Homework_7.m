%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 7
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 29/04/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 18;

%% VARIAVEIS
m = 10;             % Massa a ser isolada [kg]
fn = 5;             % Freq. de Isolamento (natural) [Hz]
z = [0.05 0.005];   % Zeta [adimensional]
Fs = 500;           % Freq. de Amostragem [Hz]
T = 200;            % Periodo [s]
dt = 1/Fs;          % Incremento de Tempo [s]
df = 1/T;           % Incremento de Frequencia [Hz]
t = 0:dt:T;         % Vetor de Tempo [s]
f = 0:df:Fs;        % Vetor de Frequencia [Hz]
w = 2*pi*f;         % Velocidade Angular [rad/s]

%% PARAMETROS DO SISTEMA
wn = 2*pi*fn;           % Freq. Natural [rad/s]
k = wn^2*m;             % Rigidez [N/m]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% EXCITACAO RANDOMICA x(t)
xt = randn(1,length(t));

%% CALCULO DE y(t) USANDO H(jw)
Xjw = fft(xt)*dt;
Hjw = [];
for st=1:2
    Hjw(st,:) = (k+1i*w*c(st))./(k-w.^2*m+1i*w*c(st));
    Yjw(st,:) = Hjw(st,:).*Xjw;
    yt(st,:) = ifft(Yjw(st,:))*Fs;
end

%% ESTIMANDO H(jw) USANDO x(t) E y(t)
N = round(length(t));
for st=1:2
    [HjwEst(st,:), fEst(st,:)] = tfestimate(xt,yt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end

%% PLOTANDO E COMPARANDO AS FRF's
% ABS(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:2
    subplot(1,2,st)
    loglog(f,abs(Hjw(st,:)),'r','linewidth', 2), hold on
    loglog(fEst(st,:),abs(HjwEst(st,:)),'--k','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([1e-1 Fs/2])
    set(gca,'fontsize',txtsize,'Xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3],'Ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% ANGULO(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:2
    subplot(1,2,st)
    semilogx(f,rad2deg(angle(Hjw(st,:))),'r','linewidth', 2), hold on
    semilogx(fEst(st,:),rad2deg(angle(HjwEst(st,:))),'--k','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([1e-1 Fs/2]), ylim([-180 5])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% OBTENCAO DE k E zeta PELA ANÁLISE DE H(jw) ESTIMADO
for st=1:2
   [HjwEstMax(st), HjwEstMaxIdx(st)] = max(abs(HjwEst(st,:)));
   zEst(st) = 1/(2*HjwEstMax(st));
   kEst(st) = m*w(HjwEstMaxIdx(st))^2;
end

%% COMPARACAO ENTRE OS RESULTADOS OBTIDOS
for st=1:2
    erroZeta(st) = abs(zEst(st)- z(st))./z(st)*100;
    errok(st) = abs(kEst(st)- k)./k*100;
    if st==1
        fprintf('Erro Zeta1: %f %%\n', erroZeta(st))
        fprintf('Erro k1: %f %%\n', errok(st))
        fprintf('\n')
    else
        fprintf('Erro Zeta2: %f %%\n', erroZeta(st))
        fprintf('Erro k2: %f %%\n', errok(st))
        fprintf('\n')
    end
end