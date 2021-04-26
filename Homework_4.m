%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK #4
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 01/04/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% VARIAVEIS
m = 1; %[kg]
k = 1e4; %[N/m]
z = [0.1 0.01 0.001]; % zeta [adimensional]
Fs = 5000; % Frequencia de amostragem [Hz]
dt = 1/Fs; %[s]
tmax = 60; % Tempo maximo avaliado [s]
t = 0:dt:tmax; %[s]

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m); %[rad/s]
fn = wn/(2*pi); %[Hz]
wd = wn*sqrt(1-z.^2); %[rad/s]
c = z*2*sqrt(k*m); %[N.s/m]

%% INPUT: f(t) ==> FUNCAO IMPULSO
fImp = zeros(1, length(t)); % Criando vetor para alocacao
fImp(t==dt) = 1; % Impulso unitario em dt
AreaImp = trapz(fImp)*dt; % Area sob curva de Impulso

% VISUALIZACAO ENTRADA IMPULSO
figure
plot(t,fImp,'b','linewidth', 3), hold on
grid on, grid minor
axis([0 3*dt 0 1.1])
xlabel('$t$ [s]')
ylabel('$f(t)$ [N]')
set(gca,'fontsize',18,'YTick',[0:0.1:1.1]) %#ok<NBRAK>

% IRF ANALITICO
htImp = []; % Criando os vetores
for st=1:3
    htImp(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t);
end

% CALCULO DA CONVOLUCAO - MÉTODO 1
xImpConv = [];
for st=1:3
	xImpConv_aux = conv(htImp(st,:),fImp)*dt;
    xImpConv = [xImpConv; xImpConv_aux(1:length(fImp))]; %#ok<*AGROW>
end

% DOMINIO DA FREQUENCIA - MÉTODO 2
FjwImp = fft(fImp)*dt;
for st=1:3
    HjwImp(st,:) = fft(htImp(st,:))*dt; %#ok<*SAGROW>
    XjwImp(st,:) = HjwImp(st,:).*FjwImp;
    xImpFreq(st,:) = ifft(XjwImp(st,:))*Fs;
end

figure % COMPARACAO: IRF ANALITICO PELA CONVOLUCAO E PELA FT
for st=1:3
    subplot(2,3,st)
    if st==1
       plot(t,htImp(st,:),'m','linewidth', 2), hold on
       axis([0 0.6 -1e-2 1e-2])
    elseif st==2
        plot(t,htImp(st,:),'k','linewidth', 2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,htImp(st,:),'b','linewidth', 2), hold on
        axis([0 60 -1e-2 1e-2])
    end
    plot(t,xImpConv(st,:)/AreaImp,'--r','linewidth', 1.8), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'IRF','Convolution Method'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18)
end
for st=1:3
    subplot(2,3,st+3)
    if st==1
       plot(t,htImp(st,:),'m','linewidth', 2), hold on
       axis([0 0.6 -1e-2 1e-2])
    elseif st==2
        plot(t,htImp(st,:),'k','linewidth', 2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,htImp(st,:),'b','linewidth', 2), hold on
        axis([0 60 -1e-2 1e-2])
    end
    plot(t,xImpFreq(st,:)/AreaImp,'--r','linewidth', 1.8), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'IRF','FT Method'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18)
end

%% INPUT: f(t) ==> FUNCAO IMPULSO SENO
Tn = 2*pi/wn; % Periodo Natural
Tp = [Tn/10 Tn 10*Tn]; % Diferentes Periodos do Impulso Seno
wt = 2*pi./Tp; % Construindo o sinal (w)
for TpAux = 1:3
    tt = 0:dt:Tp(TpAux)/2; % Construindo o sinal (t)
    xImpSeno = sin(wt(TpAux).*tt); % Construindo o sinal sin(wt)
    xImpSeno0 = zeros(1,length(t)-length(xImpSeno)); % Zerando o restante
    fSeno(TpAux,1:length(xImpSeno)) = xImpSeno; % Criando vetor para alocacao
    fSeno(TpAux,(1+length(xImpSeno)):length(t)) = xImpSeno0;
    AreaSeno(TpAux,:) = trapz(fSeno(TpAux,:))*dt;
end

% VISUALIZACAO DAS ENTRADAS IMPULSO SENO
figure
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(t,fSeno(st,:),'m','linewidth', 3), hold on
        axis([0 0.006 0 1.1])
        set(gca,'fontsize',18,'XTick',[0:0.001:0.006],'YTick',[0:0.1:1.1]) %#ok<NBRAK>
    elseif st==2
        plot(t,fSeno(st,:),'k','linewidth', 3), hold on
        axis([0 0.06 0 1.1])
        set(gca,'fontsize',18,'XTick',[0:0.01:0.06],'YTick',[0:0.1:1.1]) %#ok<NBRAK>
    else
        plot(t,fSeno(st,:),'b','linewidth', 3), hold on
        axis([0 0.6 0 1.1])
        set(gca,'fontsize',18,'XTick',[0:0.1:0.6],'YTick',[0:0.1:1.1]) %#ok<NBRAK>
    end
    title(['$T_p = ', num2str(Tp(st),'%.3f'), ' \,\textrm{[s]}$'],'FontWeight','normal')
    grid on, grid minor
    xlabel('$t$ [s]')
    ylabel('$f(t)$ [N]')
end

% CALCULO DA CONVOLUCAO - MÉTODO 1
xSenoConv = [];
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(1,:))*dt;%Tp1
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(2,:))*dt;%Tp2
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(3,:))*dt;%Tp3
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end

% DOMINIO DA FREQUENCIA - MÉTODO 2
xSenoFreq = [];
FjwSeno = fft(fSeno(1,:))*dt; %Tp1
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end
FjwSeno = fft(fSeno(2,:))*dt; %Tp2
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end
FjwSeno = fft(fSeno(3,:))*dt; %Tp3
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end

figure % COMPARACAO Tp1
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st,:),'m','linewidth', 2), hold on
    plot(t,xSenoFreq(st,:),'--k','linewidth', 1.8), hold on
    if st==1
        axis([0 0.6 -2e-5 2e-5])
    elseif st==2
        axis([0 6 -2e-5 2e-5])
    else
        axis([0 60 -2e-5 2e-5])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_1}$','FT Method, $T_{p_1}$'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18)
end
figure % COMPARACAO Tp2
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st+3,:),'k','linewidth', 2), hold on
    plot(t,xSenoFreq(st+3,:),'--r','linewidth', 1.8), hold on
    if st==1
        axis([0 0.6 -1.6e-4 1.6e-4])
    elseif st==2
        axis([0 6 -1.6e-4 1.6e-4])
    else
        axis([0 60 -1.6e-4 1.6e-4])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_2}$','FT Method, $T_{p_2}$'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18)
end
figure % COMPARACAO Tp3
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st+6,:),'b','linewidth', 2), hold on
    plot(t,xSenoFreq(st+6,:),'--k','linewidth', 1.8), hold on
    if st==1
        axis([0 0.6 -2.1e-5 1.2e-4])
    elseif st==2
        axis([0 6 -2.1e-5 1.2e-4])
    else
        axis([0 60 -2.1e-5 1.2e-4])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_3}$','FT Method, $T_{p_3}$'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18)
end