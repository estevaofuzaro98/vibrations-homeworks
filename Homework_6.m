%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 6
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 22/04/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 17;

%% VARIAVEIS
m = 1;                  % Massa [kg]
k = 1e4;                % Rigidez [N/m]
z = [0.001 0.01 0.1];   % Zeta [adimensional]
Fs = 1000;              % Freq. de Amostragem [Hz]
dt = 1/Fs;              % Incremento de Tempo [s]
T = 100;                % Periodo [s]
df = 1/T;               % Incremento de Freq. [Hz]
t = 0:dt:T;             % Vetor de Tempo [s]
f = 0:df:(Fs/2);        % Frequencia [Hz]
w = 2*pi*f;             % Freq. Angular [rad/s]
SNR = 10;               % Signal to Noise Ratio

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m);         % Freq. Natural [rad/s]
fn = wn/(2*pi);         % Freq. Natural [Hz]
wd = wn*sqrt(1-z.^2);   % Freq. Nat. Amortecida [rad/s]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% IRF E FRF ANALITICAS
ht = []; Hjw = []; % Pré-alocando os vetores
for st=1:3
    ht(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t); %#ok<*SAGROW>
    Hjw(st,:) = 1./(k-w.^2*m+1i*w*c(st));
end

%% EXCITACAO RANDOMICA f(t) PURA
ft = randn(1,length(t));
% CALCULO DE x(t) PELA CONVOLUCAO
xt = [];
for st=1:3
    x_aux = conv(ht(st,:),ft)*dt;
    xt = [xt; x_aux(1:length(ft))]; %#ok<*AGROW>
end
% ESTIMANDO H(jw) E PLOTANDO AS FRF's
N = round(length(t));
for st=1:3
    [HjwEst(st,:), fEst(st,:)] = tfestimate(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure % ABS(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
    elseif st==2
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
    else
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
    end
    loglog(fEst(st,:),abs(HjwEst(st,:)),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([1e0 1e2])
    set(gca,'fontsize',txtsize,'Xtick',[1e0 1e1 1e2 1e3],'Ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'XColor','k','YColor','k','ZColor','k','GridColor','k')
    if st==1
        axes('Position',[.16 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([0.01 0.06])
        set(gca,'Xtick',[],'Ytick',[])
    elseif st==2
        axes('Position',[.44 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([4.5 5.2]*1e-3)
        set(gca,'Xtick',[],'Ytick',[])
    else
        axes('Position',[.72 .34 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([14.6 16.8]); ylim([4.2 5.6]*1e-4)
        set(gca,'Xtick',[],'Ytick',[])
    end
end
figure % ANGULO(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'m','linewidth', 2), hold on
        xlim([1e0 25])
    elseif st==2
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'k','linewidth', 2), hold on
        xlim([1e0 50])
    else
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'b','linewidth', 2), hold on
        xlim([1e0 1e2])
    end
    semilogx(fEst(st,:),rad2deg(unwrap(angle(HjwEst(st,:)))),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([-180 5])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% FUNCAO DE COERENCIA
for st=1:3
    coh(st,:) = mscohere(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(fEst(st,:),coh(st,:),'m','linewidth', 2), hold on
    elseif st==2
        semilogx(fEst(st,:),coh(st,:),'k','linewidth', 2), hold on
    else
        semilogx(fEst(st,:),coh(st,:),'b','linewidth', 2), hold on
    end
    xlabel('$f$ [Hz]')
    ylabel('$\gamma^2$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([0 1])
    set(gca,'fontsize',txtsize,'Xtick',[1e-1 1e0 1e1 1e2 1e3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% EXCITACAO RANDOMICA f(t) COM RUIDO
ft = randn(1,length(t));
ft = ft + sqrt(1/(10^(SNR/10))*var(ft))*randn(1,length(t));
% CALCULO DE x(t) PELA CONVOLUCAO
xt = [];
for st=1:3
    x_aux = conv(ht(st,:),ft)*dt;
    xt = [xt; x_aux(1:length(ft))]; %#ok<*AGROW>
end
% ESTIMANDO H(jw) E PLOTANDO AS FRF's
N = round(length(t));
for st=1:3
    [HjwEst(st,:), fEst(st,:)] = tfestimate(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure % ABS(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
    elseif st==2
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
    else
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
    end
    loglog(fEst(st,:),abs(HjwEst(st,:)),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([1e0 1e2])
    set(gca,'fontsize',txtsize,'Xtick',[1e0 1e1 1e2 1e3],'Ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'XColor','k','YColor','k','ZColor','k','GridColor','k')
    if st==1
        axes('Position',[.16 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([0.01 0.06])
        set(gca,'Xtick',[],'Ytick',[])
    elseif st==2
        axes('Position',[.44 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([4.5 5.2]*1e-3)
        set(gca,'Xtick',[],'Ytick',[])
    else
        axes('Position',[.72 .34 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([14.6 16.8]); ylim([4.2 5.6]*1e-4)
        set(gca,'Xtick',[],'Ytick',[])
    end
end
figure % ANGULO(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'m','linewidth', 2), hold on
        xlim([1e0 25])
    elseif st==2
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'k','linewidth', 2), hold on
        xlim([1e0 50])
    else
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'b','linewidth', 2), hold on
        xlim([1e0 1e2])
    end
    semilogx(fEst(st,:),rad2deg(unwrap(angle(HjwEst(st,:)))),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([-180 5])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% FUNCAO DE COERENCIA
for st=1:3
    coh(st,:) = mscohere(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(fEst(st,:),coh(st,:),'m','linewidth', 2), hold on
    elseif st==2
        semilogx(fEst(st,:),coh(st,:),'k','linewidth', 2), hold on
    else
        semilogx(fEst(st,:),coh(st,:),'b','linewidth', 2), hold on
    end
    xlabel('$f$ [Hz]')
    ylabel('$\gamma^2$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([0 1])
    set(gca,'fontsize',txtsize,'Xtick',[1e-1 1e0 1e1 1e2 1e3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% EXCITACAO RANDOMICA f(t) COM RUIDO EM x(t)
ft = randn(1,length(t));
% CALCULO DE x(t) PELA CONVOLUCAO
xtAux = [];
for st=1:3
    x_aux = conv(ht(st,:),ft)*dt;
    xtAux = [xtAux; x_aux(1:length(ft))]; %#ok<*AGROW>
end
for st=1:3
    xt(st,:) = xtAux(st,1:length(t)) + sqrt(1/(10^(SNR/10))*var(xtAux(st,:)))*randn(1,length(t));
end
% ESTIMANDO H(jw) E PLOTANDO AS FRF's
N = round(length(t));
for st=1:3
    [HjwEst(st,:), fEst(st,:)] = tfestimate(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure % ABS(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
    elseif st==2
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
    else
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
    end
    loglog(fEst(st,:),abs(HjwEst(st,:)),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([1e0 1e2])
    set(gca,'fontsize',txtsize,'Xtick',[1e0 1e1 1e2 1e3],'Ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2],'XColor','k','YColor','k','ZColor','k','GridColor','k')
    if st==1
        axes('Position',[.16 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([0.01 0.06])
        set(gca,'Xtick',[],'Ytick',[])
    elseif st==2
        axes('Position',[.44 .61 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([15.8 16]); ylim([4.5 5.2]*1e-3)
        set(gca,'Xtick',[],'Ytick',[])
    else
        axes('Position',[.72 .34 .062 .27]); box on
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 3), hold on
        loglog(fEst(st,:),abs(HjwEst(st,:)),'r--','linewidth', 3), hold on
        xlim([14.6 16.8]); ylim([4.2 5.6]*1e-4)
        set(gca,'Xtick',[],'Ytick',[])
    end
end
figure % ANGULO(X/F) PELA FREQUENCIA
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'m','linewidth', 2), hold on
        xlim([1e0 25])
    elseif st==2
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'k','linewidth', 2), hold on
        xlim([1e0 50])
    else
        semilogx(f,rad2deg(unwrap(angle(Hjw(st,:)))),'b','linewidth', 2), hold on
        xlim([1e0 1e2])
    end
    semilogx(fEst(st,:),rad2deg(unwrap(angle(HjwEst(st,:)))),'--r','linewidth', 2), hold on
    xlabel('$f$ [Hz]')
    legend({'$H(j\omega)$ Analytic','$H(j\omega)$ Estimated'},'Location','southwest','fontsize',lgndsize)
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([-180 5])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% FUNCAO DE COERENCIA
for st=1:3
    coh(st,:) = mscohere(ft,xt(st,:),[],[],N,Fs); %#ok<*SAGROW>
end
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        semilogx(fEst(st,:),coh(st,:),'m','linewidth', 2), hold on
    elseif st==2
        semilogx(fEst(st,:),coh(st,:),'k','linewidth', 2), hold on
    else
        semilogx(fEst(st,:),coh(st,:),'b','linewidth', 2), hold on
    end
    xlabel('$f$ [Hz]')
    ylabel('$\gamma^2$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([0 1])
    set(gca,'fontsize',txtsize,'Xtick',[1e-1 1e0 1e1 1e2 1e3],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end