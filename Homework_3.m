%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 3
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 25/03/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 20;

%% VARIAVEIS
m = 1;                  % Massa [kg]
k = 1e4;                % Rigidez [N/m]
z = [0.001 0.01 0.1];   % Zeta [adimensional]
dt = 0.001;             % Incremento de Tempo [s]
df = 0.001;             % Incremento de Frequencia [Hz]
fs = 1/dt;              % Freq. de Amostragem [Hz]
t = 0:dt:200;           % Vetor de Tempo [s]
f = 0:df:(fs/2);        % Frequencia [Hz]
w = 2*pi*f;             % Freq. Angular [rad/s]

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m);         % Freq. Natural [rad/s]
fn = wn/(2*pi);         % Freq. Natural [Hz]
wd = wn*sqrt(1-z.^2);   % Freq. Nat. Amortecida [rad/s]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% IRF [h(t)]
ht = []; Hjw = []; env = []; % Pré-alocando os vetores
for st=1:3
    ht(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t); %#ok<*SAGROW>
    Hjw(st,:) = 1./(k-w.^2*m+1i*w*c(st));
    env(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t);
end

%% TRANSFORMADA DE FOURIER DE h(t)
Ts = [160 40 10];   % T Amostragem para cada zeta [s]
DiscFreq = 1./Ts;   % Freq. Discretizada [Hz]

% ABS(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        Fs = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    elseif st==2
        Fs = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    else
        Fs = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    end
    HjwFT = fft(ht(st,1:length(Fs)))*dt;
    loglog(Fs(1:fs/(2*DiscFreq(st))),abs(HjwFT(1:fs/(2*DiscFreq(st)))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$\left|H(j\omega)\right|$','$\left|FT[h(t)]\right|$'},'Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'Xtick',[1e-1 1e0 1e1 1e2 1e3],'Ytick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% ANGULO(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        Fs = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'m','linewidth', 2), hold on
    elseif st==2
        Fs = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'k','linewidth', 2), hold on
    else
        Fs = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'b','linewidth', 2), hold on
    end
    HjwFT = fft(ht(st,1:length(Fs)))*dt;
    plot(Fs(1:fs/(2*DiscFreq(st))),rad2deg(angle(HjwFT(1:fs/(2*DiscFreq(st))))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    legend({'$\angle H(j\omega)$','$\angle FT[h(t)]$'},'Location','northeast','fontsize',lgndsize)
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([0 80])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% DOUBLE-SIDED SPECTRUM [DSS] PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    DSfreq = 0:df:fs;
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DSfreq));
    subplot(1,3,st)
    if st==1
        Fs = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'m','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    elseif st==2
        Fs = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'k','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    else
        Fs = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'b','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    end
    HjwFT = fft(ht(st,1:length(Fs)))*dt;
    semilogy(Fs,abs(HjwFT),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'DSS of H(j$\omega$)','DSS of FT[$h(t)$]'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',txtsize,'Ytick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% ANGULO DO DSS PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    DSfreq = 0:df:fs;
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DSfreq));
    subplot(1,3,st)
    if st==1
        Fs = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'m','linewidth', 2), hold on
    elseif st==2
        Fs = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'k','linewidth', 2), hold on
    else
        Fs = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'b','linewidth', 2), hold on
    end
    HjwFT = fft(ht(st,1:length(Fs)))*dt;
    plot(Fs,rad2deg(unwrap(angle(HjwFT))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    legend({'DSS of H(j$\omega$)','DSS of FT[$h(t]$)'},'Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([-360 0])
    set(gca,'fontsize',txtsize,'Ytick',[-360 -270 -180 -90 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
%% PLOT DE h(t) E DA IFT[H(jw)] PELO TEMPO
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
    DStime = 0:1/fs:1/df;
for st=1:3
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DStime));
    htFT = ifft(DSSHjw)*fs;
    subplot(1,3,st)
    if st==1
        plot(t,ht(st,:),'m','linewidth', 2), hold on
        plot(DStime,htFT,'--r','linewidth', 2.2), hold on
        axis([0 60 -1e-2 1e-2])
    elseif st == 2
        plot(t,ht(st,:),'k','linewidth', 2), hold on
        plot(DStime,htFT,'--r','linewidth', 2.2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,ht(st,:),'b','linewidth', 2), hold on
        plot(DStime,htFT,'--r','linewidth', 2.2), hold on
        axis([0 0.6 -1e-2 1e-2])
    end
    xlabel('$t$ [s]')
    ylabel('$h(t)$ [m/N.s]')
    legend({'$h(t)$','IFT[H(j$\omega$)]'},'Location','southeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    grid on; grid minor;
end
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DStime));
    htFT = ifft(DSSHjw)*fs;
    subplot(1,3,st)
    if st==1
        plot(t,ht(st,:),'m','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.2 -1e-2 1e-2])
    elseif st == 2
        plot(t,ht(st,:),'k','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.2 -1e-2 1e-2])
    else
        plot(t,ht(st,:),'b','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.2 -1e-2 1e-2])
    end
    xlabel('$t$ [s]')
    ylabel('$h(t)$ [m/N.s]')
    legend({'$h(t)$','IFT[H(j$\omega$)]'},'Location','southeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    grid on; grid minor;
end