%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 8
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 06/05/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 18;

%% VARIAVEIS
ms = 1;         % Massa principal [kg]
ks = 1e4;       % Rigidez principal [N/m]
zs = 0.001;     % Zeta principal [adimensional]
mu = 0.1;       % mu = ma/ms [adimensional]
Fs = 1000;      % Freq. de Amostragem [Hz]
T = 60;         % Periodo [s]
dt = 1/Fs;      % Incremento de Tempo [s]
df = 1/T;       % Incremento de Frequencia [Hz]
t = 0:dt:T;     % Vetor de Tempo [s]
f = 0:df:Fs;    % Vetor de Frequencia [Hz]
w = 2*pi*f;     % Velocidade Angular [rad/s]

%% PARAMETROS DO SISTEMA
wn = sqrt(ks/ms);               % Freq. Natural [rad/s]
cs = 2*zs*sqrt(ks*ms);          % Amortecimento principal [N.s/m]
wa = wn/(1+mu);                 % Freq. Nat. do absorvedor [Hz]
ma = ms*mu;                     % Massa do absorvedor [kg]
ka = ma*wa^2;                   % Rigidez do absorvedor [N/m]
za = sqrt((3/8)*(mu/(1+mu)^3)); % Zeta do absorvedor [adimensional]
ca = 2*za*sqrt(ka*ma);          % Amortecimento do absorvedor [N.s/m]

%% CALCULO DAS FRF's Hs E Ha
M = [ms 0; 0 ma];
K = [ks+ka -ka; -ka ka];    % Modelo
C = [cs+ca -ca; -ca ca];    % matricial
F = [1; 0];
st=0;
for fAux=0:df:Fs
   st=st+1;
   wAux = 2*pi*fAux;
   D = K - wAux.^2*M + 1i*wAux*C;
   H = D\F;
   Hs(st) = H(1);  %#ok<*SAGROW>
   Ha(st) = H(2);
end

% PLOTANDO Hs E Ha
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
loglog(f,abs(Hs),'r','linewidth', 2), hold on
loglog(f,abs(Ha),'k--','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance [m/N]')
legend({'$H_s$','$H_a$'},'Location','northeast','fontsize',lgndsize)
grid on, grid minor
axis([1e-1 Fs/2 1e-7 2e-3])
set(gca,'fontsize',txtsize,'Xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3],'Ytick',[1-7 1e-6 1e-5 1e-4 1e-3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
axes('Position',[.25 .25 .3 .4]); box on
loglog(f,abs(Hs),'r','linewidth',2), hold on, grid on
loglog(f,abs(Ha),'k--','linewidth', 2), hold on
xlim([5 25]); ylim([1e-4 2e-3]);
grid on, grid minor
set(gca,'fontsize',txtsize-8,'Xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3],'Ytick',[1-7 1e-6 1e-5 1e-4 1e-3],'XColor','k','YColor','k','ZColor','k','GridColor','k')

figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
semilogx(f,rad2deg(unwrap(angle(Hs))),'r','linewidth', 2), hold on
semilogx(f,rad2deg(unwrap(angle(Ha))),'k--','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance $\phi$ [$^{\circ}$]')
legend({'$H_s$','$H_a$'},'Location','northeast','fontsize',lgndsize)
grid on, grid minor
axis([1e-1 Fs/2 -315 5])
set(gca,'fontsize',txtsize,'Xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
axes('Position',[.2 .3 .3 .4]); box on
semilogx(f,rad2deg(unwrap(angle(Hs))),'r','linewidth',2), hold on, grid on
semilogx(f,rad2deg(unwrap(angle(Ha))),'k--','linewidth', 2), hold on
xlim([5 25]); ylim([-315 5]);
grid on, grid minor
set(gca,'fontsize',txtsize-8,'Xtick',[1e-2 1e-1 1e0 1e1 1e2 1e3],'XColor','k','YColor','k','ZColor','k','GridColor','k')

%% EXCITACAO RANDOMICA x(t)
x = randn(1,length(t));

%% CALCULO DE y(t) USANDO H(jw)
X = fft(x)*dt;
Xs = Hs.*X;
Xa = Ha.*X;
xs = ifft(Xs)*Fs;
xa = ifft(Xa)*Fs;

% PLOTANDO xs, xa
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,xa,'r','linewidth', 2), hold on
plot(t,xs,'k','linewidth', 2), hold on
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$x(t)$ [m]','Interpreter','Latex')
legend({'$x_s(t)$','$x_a(t)$'},'location','northeast','fontsize',lgndsize)
xlim([0 50]); ylim([-2.5e-4 2.5e-4])
grid on, grid minor
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')

%% ESTIMANDO Ha e Hs
N = round(length(t));
[HsEst, fsEst] = tfestimate(x,xs,[],[],N,Fs);
[HaEst, faEst] = tfestimate(x,xa,[],[],N,Fs);

%% PLOTANDO E COMPARANDO AS FRF's
% ABS(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
subplot(1,2,1) %Hs
loglog(f,abs(Hs),'r','linewidth', 2), hold on
loglog(fsEst,abs(HsEst),'--k','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance [m/N]')
legend({'$H_s(j\omega)$ Analytic','$H_s(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
grid on, grid minor
xlim([1e-1 Fs/2])
set(gca,'fontsize',txtsize,'Ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
subplot(1,2,2) %Ha
loglog(f,abs(Ha),'r','linewidth', 2), hold on
loglog(faEst,abs(HaEst),'--k','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance [m/N]')
legend({'$H_a(j\omega)$ Analytic','$H_a(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
grid on, grid minor
xlim([1e-1 Fs/2])
set(gca,'fontsize',txtsize,'Ytick',[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3],'XColor','k','YColor','k','ZColor','k','GridColor','k')

% ANGULO(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
subplot(1,2,1) %Hs
semilogx(f,rad2deg(unwrap(angle(Hs))),'r','linewidth', 2), hold on
semilogx(fsEst,rad2deg(unwrap(angle(HsEst))),'--k','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance $\phi$ [$^{\circ}$]')
legend({'$H_s(j\omega)$ Analytic','$H_s(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
grid on, grid minor
axis([1e-1 1e2 -200 10])
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
subplot(1,2,2) %Ha
semilogx(f,rad2deg(unwrap(angle(Ha))),'r','linewidth', 2), hold on
semilogx(faEst,rad2deg(unwrap(angle(HaEst))),'--k','linewidth', 2), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance $\phi$ [$^{\circ}$]')
legend({'$H_a(j\omega)$ Analytic','$H_a(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
grid on, grid minor
axis([1e-1 80 -350 10])
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')