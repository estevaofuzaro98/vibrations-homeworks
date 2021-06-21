%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 4
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 01/04/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 18;

%% VARIAVEIS
m = 1;                  % Massa [kg]
k = 1e4;                % Rigidez [N/m]
z = [0.1 0.01 0.001];   % Zeta [adimensional]
Fs = 5000;              % Freq. de Amostragem [Hz]
dt = 1/Fs;              % Incremento de Tempo [s]
t = 0:dt:60;            % Vetor de Tempo [s]

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m);         % Freq. Natural [rad/s]
fn = wn/(2*pi);         % Freq. Natural [Hz]
wd = wn*sqrt(1-z.^2);   % Freq. Nat. Amortecida [rad/s]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% INPUT: f(t) ==> FUNCAO IMPULSO
fImp = zeros(1, length(t)); % Criando vetor para alocacao
fImp(t==dt) = 1;            % Impulso unitario em dt
AreaImp = trapz(fImp)*dt;   % Area sob curva de Impulso

% VISUALIZACAO ENTRADA IMPULSO
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,fImp,'b','linewidth', 3), hold on
grid on, grid minor
axis([0 3*dt 0 1.1])
xlabel('$t$ [s]')
ylabel('$f(t)$ [N]')
set(gca,'fontsize',txtsize,'YTick',[0:0.2:1.1],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>

% IRF ANALITICO
htImp = []; % Criando os vetores
for st=1:3
    htImp(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t);
end

% CALCULO DA CONVOLUCAO - METODO 1
xImpConv = [];
for st=1:3
    xImpConv_aux = conv(htImp(st,:),fImp)*dt;
    xImpConv = [xImpConv; xImpConv_aux(1:length(fImp))]; %#ok<*AGROW>
end

% DOMINIO DA FREQUENCIA - METODO 2
FjwImp = fft(fImp)*dt;
for st=1:3
    HjwImp(st,:) = fft(htImp(st,:))*dt; %#ok<*SAGROW>
    XjwImp(st,:) = HjwImp(st,:).*FjwImp;
    xImpFreq(st,:) = ifft(XjwImp(st,:))*Fs;
end

figure % COMPARACAO: IRF ANALITICO PELA CONVOLUCAO E FT
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
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
    plot(t,xImpConv(st,:)/AreaImp,'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'IRF','Convolution Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
figure % COMPARACAO: IRF ANALITICO PELA CONVOLUCAO E FT
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
       plot(t,htImp(st,:),'m','linewidth', 2), hold on
       axis([0 0.6 -1e-2 1e-2])
    elseif st==2
        plot(t,htImp(st,:),'k','linewidth', 2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,htImp(st,:),'b','linewidth', 3), hold on
        axis([0 60 -1e-2 1e-2])
    end
    plot(t,xImpFreq(st,:)/AreaImp,'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'IRF','FT Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

figure % COMPARACAO: IRF PELA CONVOLUCAO E FT
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
       plot(t,xImpConv(st,:)/AreaImp,'m','linewidth', 2), hold on
       axis([0 0.6 -1e-2 1e-2])
    elseif st==2
        plot(t,xImpConv(st,:)/AreaImp,'k','linewidth', 2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,xImpConv(st,:)/AreaImp,'b','linewidth', 2), hold on
        axis([0 60 -1e-2 1e-2])
    end
    plot(t,xImpFreq(st,:)/AreaImp,'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method','FT Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% INPUT: f(t) ==> FUNCAO IMPULSO SENO
Tn = 2*pi/wn;           % Periodo Natural
Tp = [Tn/10 Tn 10*Tn];  % Diferentes Periodos de Impulso Seno
wt = 2*pi./Tp;          % Construindo o sinal (w)
for TpAux = 1:3
    tt = 0:dt:Tp(TpAux)/2;                              % Construindo o sinal (t)
    xImpSeno = sin(wt(TpAux).*tt);                      % Construindo o sinal sin(wt)
    xImpSeno0 = zeros(1,length(t)-length(xImpSeno));    % Zerando o restante
    fSeno(TpAux,1:length(xImpSeno)) = xImpSeno;
    fSeno(TpAux,(1+length(xImpSeno)):length(t)) = xImpSeno0;
    AreaSeno(TpAux,:) = trapz(fSeno(TpAux,:))*dt;
end

% VISUALIZACAO DAS ENTRADAS IMPULSO SENO
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(t,fSeno(st,:),'m','linewidth', 3), hold on
        axis([0 0.006 0 1.1])
        set(gca,'fontsize',txtsize,'XTick',[0:0.002:0.006],'YTick',[0:0.2:1.1],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>
    elseif st==2
        plot(t,fSeno(st,:),'k','linewidth', 3), hold on
        axis([0 0.06 0 1.1])
        set(gca,'fontsize',txtsize,'XTick',[0:0.02:0.06],'YTick',[0:0.2:1.1],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>
    else
        plot(t,fSeno(st,:),'b','linewidth', 3), hold on
        axis([0 0.6 0 1.1])
        set(gca,'fontsize',txtsize,'XTick',[0:0.2:0.6],'YTick',[0:0.2:1.1],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>
    end
    title(['$T_p = ', num2str(Tp(st),'%.3f'), ' \,\textrm{[s]}$'],'FontWeight','normal')
    grid on, grid minor
    xlabel('$t$ [s]')
    ylabel('$f(t)$ [N]')
end

% CALCULO DA CONVOLUCAO - METODO 1
xSenoConv = [];
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(1,:))*dt;    % Tp1
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(2,:))*dt;    % Tp2
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end
for st=1:3 
    xSenoConv_aux = conv(htImp(st,:),fSeno(3,:))*dt;    % Tp3
    xSenoConv = [xSenoConv; xSenoConv_aux(1,1:length(fSeno))];
end

% DOMINIO DA FREQUENCIA - METODO 2
xSenoFreq = [];
FjwSeno = fft(fSeno(1,:))*dt;   % Tp1
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end
FjwSeno = fft(fSeno(2,:))*dt;   % Tp2
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end
FjwSeno = fft(fSeno(3,:))*dt;   % Tp3
for st=1:3
    HjwSeno(st,:) = fft(htImp(st,:))*dt;
    XjwSeno(st,:) = HjwSeno(st,:).*FjwSeno;
    xSenoFreq_aux = ifft(XjwSeno(st,:))*Fs;
    xSenoFreq = [xSenoFreq; xSenoFreq_aux(1,:)];
end

figure % COMPARACAO Tp1
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st,:),'m','linewidth', 2), hold on
    plot(t,xSenoFreq(st,:),'--k','linewidth', 1.4), hold on
    if st==1
        axis([0 0.6 -2e-5 2e-5])
    elseif st==2
        axis([0 6 -2e-5 2e-5])
    else
        axis([0 60 -2e-5 2e-5])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_1}$','FT Method, $T_{p_1}$'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
figure % COMPARACAO Tp2
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st+3,:),'k','linewidth', 2), hold on
    plot(t,xSenoFreq(st+3,:),'--r','linewidth', 1.4), hold on
    if st==1
        axis([0 0.6 -1.6e-4 1.6e-4])
    elseif st==2
        axis([0 6 -1.6e-4 1.6e-4])
    else
        axis([0 60 -1.6e-4 1.6e-4])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_2}$','FT Method, $T_{p_2}$'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
figure % COMPARACAO Tp3
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    plot(t,xSenoConv(st+6,:),'b','linewidth', 2), hold on
    plot(t,xSenoFreq(st+6,:),'--k','linewidth', 1.4), hold on
    if st==1
        axis([0 0.6 -2.1e-5 1.2e-4])
    elseif st==2
        axis([0 6 -2.1e-5 1.2e-4])
    else
        axis([0 60 -2.1e-5 1.2e-4])
    end
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method, $T_{p_3}$','FT Method, $T_{p_3}$'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% INPUT: f(t) ==> FUNCAO RANDOMICA
fRand = randn(1,length(t));

% VISUALIZACAO ENTRADA RANDOMICA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,fRand,'b','linewidth', 3), hold on
grid on, grid minor
axis([0 60 -5 5])
xlabel('$t$ [s]')
ylabel('$f(t)$ [N]')
set(gca,'fontsize',txtsize,'YTick',[-6:2:6],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>

% CALCULO DA CONVOLUCAO - METODO 1
xRandConv = [];
for st=1:3
    xRandConv_aux = conv(htImp(st,:),fRand)*dt;
    xRandConv = [xRandConv; xRandConv_aux(1:length(fRand))]; %#ok<*AGROW>
end

% DOMINIO DA FREQUENCIA - METODO 2
FjwRand = fft(fRand)*dt;
for st=1:3
    HjwRand(st,:) = fft(htImp(st,:))*dt; %#ok<*SAGROW>
    XjwRand(st,:) = HjwRand(st,:).*FjwRand;
    xRandFreq(st,:) = ifft(XjwRand(st,:))*Fs;
end

figure % COMPARACAO: IRF ANALITICO PELA CONVOLUCAO E FT
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
       plot(t,xRandConv(st,:),'m','linewidth', 2), hold on
       axis([0 60 -1e-4 1e-4])
    elseif st==2
        plot(t,xRandConv(st,:),'k','linewidth', 2), hold on
        axis([0 60 -3e-4 3e-4])
    else
        plot(t,xRandConv(st,:),'b','linewidth', 2), hold on
        axis([0 60 -7e-4 7e-4])
    end
    plot(t,xRandFreq(st,:),'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method','FT Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% INPUT: f(t) ==> FUNCAO CHIRP [1 a 100 Hz]
fChirp1s = chirp(t,1,1,100);    % Chirp de 1 segundo
fChirp10s = chirp(t,1,10,100);  % Chirp de 10 segundos

% VISUALIZACAO ENTRADAS CHIRP
figure  % 1 Segundo
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,fChirp1s,'b','linewidth', 3), hold on
grid on, grid minor
axis([0 1 -1 1])
xlabel('$t$ [s]')
ylabel('$f(t)$ [N]')
set(gca,'fontsize',txtsize,'YTick',[-5:1:5],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>
figure  % 10 Segundos
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,fChirp10s,'b','linewidth', 3), hold on
grid on, grid minor
axis([0 10 -1 1])
xlabel('$t$ [s]')
ylabel('$f(t)$ [N]')
set(gca,'fontsize',txtsize,'YTick',[-5:1:5],'XColor','k','YColor','k','ZColor','k','GridColor','k') %#ok<NBRAK>

% CALCULO DA CONVOLUCAO CHIRP 1s - METODO 1
xChirp1sConv = [];
for st=1:3
    xChirp1sConv_aux = conv(htImp(st,:),fChirp1s)*dt;
    xChirp1sConv = [xChirp1sConv; xChirp1sConv_aux(1:length(fChirp1s))]; %#ok<*AGROW>
end
% CALCULO DA CONVOLUCAO CHIRP 10s - METODO 1
xChirp10sConv = [];
for st=1:3
    xChirp10sConv_aux = conv(htImp(st,:),fChirp10s)*dt;
    xChirp10sConv = [xChirp10sConv; xChirp10sConv_aux(1:length(fChirp10s))]; %#ok<*AGROW>
end

% DOMINIO DA FREQUENCIA CHIRP 1s - METODO 2
FjwChirp1s = fft(fChirp1s)*dt;
for st=1:3
    HjwChirp1s(st,:) = fft(htImp(st,:))*dt; %#ok<*SAGROW>
    XjwChirp1s(st,:) = HjwChirp1s(st,:).*FjwChirp1s;
    xChirp1sFreq(st,:) = ifft(XjwChirp1s(st,:))*Fs;
end
% DOMINIO DA FREQUENCIA CHIRP 10s - METODO 2
FjwChirp10s = fft(fChirp10s)*dt;
for st=1:3
    HjwChirp10s(st,:) = fft(htImp(st,:))*dt; %#ok<*SAGROW>
    XjwChirp10s(st,:) = HjwChirp10s(st,:).*FjwChirp10s;
    xChirp10sFreq(st,:) = ifft(XjwChirp10s(st,:))*Fs;
end

figure % COMPARACAO CONVOLUCAO E FT CHIRP 1s
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
       plot(t,xChirp1sConv(st,:),'m','linewidth', 2), hold on
       axis([0 1 -3.2e-4 3.2e-4])
    elseif st==2
        plot(t,xChirp1sConv(st,:),'k','linewidth', 2), hold on
        axis([0 6 -5.5e-4 5.5e-4])
    else
        plot(t,xChirp1sConv(st,:),'b','linewidth', 2), hold on
        axis([0 50 -6.2e-4 6.2e-4])
    end
    plot(t,xChirp1sFreq(st,:),'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method','FT Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
figure % COMPARACAO CONVOLUCAO E FT CHIRP 10s
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
       plot(t,xChirp10sConv(st,:),'m','linewidth', 2), hold on
       xlim([0 10])
    elseif st==2
        plot(t,xChirp10sConv(st,:),'k','linewidth', 2), hold on
        xlim([0 10])
    else
        plot(t,xChirp10sConv(st,:),'b','linewidth', 2), hold on
        xlim([0 60])
    end
    plot(t,xChirp10sFreq(st,:),'--r','linewidth', 1.4), hold on
    xlabel('$t$ [s]')
    ylabel('$x(t)$ [m]')
    legend({'Convolution Method','FT Method'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end