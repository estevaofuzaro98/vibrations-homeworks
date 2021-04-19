%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK #3
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 25/03/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% VARIAVEIS
m = 1; %[kg]
k = 1e4; %[N/m]
z = [0.001 0.01 0.1]; %[adimensional]
dt = 0.001; %[s]
df = 0.001; %[Hz]
fs = 1/dt; %freq de amostragem [Hz]
t = 0:dt:200; %[s]
f = 0:df:(fs/2); %[Hz]
w = 2*pi*f; %[rad/s]
figc = 1; % contador de figuras

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m); %[rad/s]
fn = wn/(2*pi); %[Hz]
wd = wn*sqrt(1-z.^2); %[rad/s]
c = z*2*sqrt(k*m); %[Ns/m]

%% IRF [h(t)]
ht = []; Hjw = []; env = []; % Criando os vetores
for st=1:3
    ht = [ht; 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t)]; %#ok<*AGROW>
    Hjw = [Hjw; 1./(k-w.^2*m+1i*w*c(st))];
    env = [env; 1/(m*wd(st))*exp(-z(st)*wn*t)];
end

%% TRANSFORMADA DE FOURIER DE h(t)
% ESCOLHA DOS PERIODOS DE AMOSTRAGEM PARA CADA ZETA
T_am = [160 40 10]; %[s]
DiscFreq = 1./T_am;

% ABS(X/F) PELA FREQUENCIA
figure(figc); figc = figc + 1;
for st=1:3
    subplot(2,3,st)
    if st==1
        f_am = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    elseif st==2
        f_am = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    else
        f_am = 0:DiscFreq(st):fs;
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
        axis([1e-1 0.5e3 1e-7 1e-1])
    end
    HjwFT = fft(ht(st,1:length(f_am)))*dt;
    loglog(f_am(1:fs/(2*DiscFreq(st))),abs(HjwFT(1:fs/(2*DiscFreq(st)))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'$\left|H(j\omega)\right|$','$\left|FT[h(t)]\right|$'},'Location','southwest','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18,'Xtick',[1e-1 1e0 1e1 1e2],'Ytick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% ANGULO(X/F) PELA FREQUENCIA
for st=1:3
    subplot(2,3,st+3)
    if st==1
        f_am = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'m','linewidth', 2), hold on
    elseif st==2
        f_am = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'k','linewidth', 2), hold on
    else
        f_am = 0:DiscFreq(st):fs;
        plot(f,rad2deg(angle(Hjw(st,:))),'b','linewidth', 2), hold on
    end
    HjwFT = fft(ht(st,1:length(f_am)))*dt;
    plot(f_am(1:fs/(2*DiscFreq(st))),rad2deg(angle(HjwFT(1:fs/(2*DiscFreq(st))))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    legend({'$\angle H(j\omega)$','$\angle FT[h(t)]$'},'Location','northeast','fontsize',15)
    ylabel('$\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    xlim([0 80])
    set(gca,'fontsize',18,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% DOUBLE-SIDED SPECTRUM [DSS] PELA FREQUENCIA
figure(figc); figc = figc + 1;
for st=1:3
    DSfreq = 0:df:fs;
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DSfreq));
    subplot(2,3,st)
    if st==1
        f_am = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'m','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    elseif st==2
        f_am = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'k','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    else
        f_am = 0:DiscFreq(st):fs;
        semilogy(DSfreq,abs(DSSHjw),'b','linewidth', 2), hold on
        axis([0 fs 1e-7 1e-1])
    end
    HjwFT = fft(ht(st,1:length(f_am)))*dt;
    semilogy(f_am,abs(HjwFT),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend({'DSS of H(j$\omega$)','FT[$h(t)$]'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    set(gca,'fontsize',18,'Ytick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
% ANGULO DO DSS PELA FREQUENCIA
for st=1:3
    DSfreq = 0:df:fs;
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DSfreq));
    subplot(2,3,st+3)
    if st==1
        f_am = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'m','linewidth', 2), hold on
    elseif st==2
        f_am = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'k','linewidth', 2), hold on
    else
        f_am = 0:DiscFreq(st):fs;
        plot(DSfreq,rad2deg(unwrap(angle(DSSHjw))),'b','linewidth', 2), hold on
    end
    HjwFT = fft(ht(st,1:length(f_am)))*dt;
    plot(f_am,rad2deg(unwrap(angle(HjwFT))),'r--','linewidth', 2.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('$\phi$ [$^{\circ}$]')
    legend({'DSS of H(j$\omega$)','FT[$h(t]$)'},'Location','northeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    grid on, grid minor
    ylim([-360 0])
    set(gca,'fontsize',18,'Ytick',[-360 -270 -180 -90 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
%% PLOT DE h(t) E DA IFT[H(jw)] PELO TEMPO
figure(figc); figc = figc + 1;
    DStime = 0:1/fs:1/df;
for st=1:3
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DStime));
    htFT = ifft(DSSHjw)*fs;
    subplot(2,3,st)
    if st==1
        plot(t,ht(st,:),'m','linewidth', 1.8), hold on
        plot(DStime,htFT,'--r','linewidth', 2), hold on
        axis([0 60 -1e-2 1e-2])
    elseif st == 2
        plot(t,ht(st,:),'k','linewidth', 1.8), hold on
        plot(DStime,htFT,'--r','linewidth', 2), hold on
        axis([0 6 -1e-2 1e-2])
    else
        plot(t,ht(st,:),'b','linewidth', 1.8), hold on
        plot(DStime,htFT,'--r','linewidth', 2), hold on
        axis([0 0.6 -1e-2 1e-2])
    end
    xlabel('$t$ [s]')
    ylabel('$h(t)$ [m/N.s]')
    legend({'$h(t)$','IFT[H(j$\omega$)]'},'Location','southeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    set(gca,'fontsize',18,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    grid on; grid minor;
end
for st=1:3
    DSSHjw = [Hjw(st,:) fliplr(conj(Hjw(st,:)))];
    DSSHjw = DSSHjw(1:length(DStime));
    htFT = ifft(DSSHjw)*fs;
    subplot(2,3,st+3)
    if st==1
        plot(t,ht(st,:),'m','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.4 -1e-2 1e-2])
    elseif st == 2
        plot(t,ht(st,:),'k','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.4 -1e-2 1e-2])
    else
        plot(t,ht(st,:),'b','linewidth', 3), hold on
        plot(DStime,htFT,'--r','linewidth', 3), hold on
        axis([0 0.4 -1e-2 1e-2])
    end
    xlabel('$t$ [s]')
    ylabel('$h(t)$ [m/N.s]')
    legend({'$h(t)$','IFT[H(j$\omega$)]'},'Location','southeast','fontsize',15)
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal')
    set(gca,'fontsize',18,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    grid on; grid minor;
end