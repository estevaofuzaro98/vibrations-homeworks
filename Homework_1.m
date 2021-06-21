%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 1
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 11/03/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 20;

%% VARIAVEIS
m = 1;                  % Massa [kg]
k = 1e4;                % Rigidez [N/m]
z = [0.1 0.01 0.001];   % Zeta [adimensional]
t = 0:0.001:100;        % Vetor de Tempo [s]
f = 0:0.001:1000;       % Frequencia [Hz]
w = 2*pi*f;             % Freq. Angular [rad/s]

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m);         % Freq. Natural [rad/s]
wd = wn*sqrt(1-z.^2);   % Freq. Nat. Amortecida [rad/s]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% IRF [h(t)] COM ENVELOPE E FRF [H(jw)] - RECEPTANCIA
ht = []; Hjw = []; env = []; % Pré-alocando os vetores
for st=1:3
    ht(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t); %#ok<*SAGROW>
    Hjw(st,:) = 1./(k-w.^2*m+1i*w*c(st));
    env(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t);
end

%% h(t) COM ENVELOPE PELO TEMPO 
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(t,ht(st,:),'m','linewidth', 1.5), hold on
    elseif st==2
        plot(t,ht(st,:),'k','linewidth', 1.5), hold on
    else
        plot(t,ht(st,:),'b','linewidth', 1.5), hold on
    end
    plot(t,env(st,:),'--r','linewidth', 1.6)
    plot(t,-env(st,:),'--r','linewidth', 1.6)
    xlabel('$t$ [s]')
    ylabel('$h(t)$ [m/N.s]')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal');
    grid on; grid minor;
    legend('$h(t)$','Envelope','','fontsize',lgndsize)
    if st==1
        axis([0 1 -1e-2 1e-2])
    elseif st==2
        axis([0 10 -1e-2 1e-2])
    else
        axis([0 100 -1e-2 1e-2])
    end
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

%% RECEPTANCIA
% ABS(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Hjw(st,:)),'m','linewidth', 2), hold on
        axis([1e-1 1e2 1e-6 1e-2])
    elseif st==2
        loglog(f,abs(Hjw(st,:)),'k','linewidth', 2), hold on
        axis([1e-1 1e2 1e-6 1e-1])
    else
        loglog(f,abs(Hjw(st,:)),'b','linewidth', 2), hold on
        axis([1e-1 1e2 1e-6 1e0])
    end
    loglog(f,1/k*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Rigidez
    loglog(f,1./(w.^2*m),'--','linewidth', 1.8), hold on % Massa
    xlabel('$f$ [Hz]')
    ylabel('Receptance [m/N]')
    legend('$\left| X/F \right|$','Stiffness Line','Mass Line','Location','southwest','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ANGULO(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(f,rad2deg(angle(Hjw(st,:))),'m','linewidth', 2), hold on
    elseif st==2
        plot(f,rad2deg(angle(Hjw(st,:))),'k','linewidth', 2), hold on
    else
        plot(f,rad2deg(angle(Hjw(st,:))),'b','linewidth', 2), hold on
    end
    plot(f,rad2deg(-pi/2)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
    xlabel('$f$ [Hz]')
    ylabel('Receptance $\phi$ [$^{\circ}$]')
    legend('$\angle \left(X/F\right)$','Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), 
    grid on, grid minor, xlim([0 50]), ylim([-185 5])
    set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ABS(X/F) PELA FREQUENCIA E ANGULO(X/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
loglog(f,abs(Hjw(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
loglog(f,abs(Hjw(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
loglog(f,abs(Hjw(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
loglog(f,1/k*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Rigidez
loglog(f,1./(w.^2*m),'--','linewidth', 1.8), hold on % Massa
xlabel('$f$ [Hz]')
ylabel('Receptance [m/N]')
legend({['$\left| X/F \right|, \, \zeta = ', num2str(z(1)),' $'],['$\left| X/F \right|, \, \zeta = ', num2str(z(2)),' $'],['$\left| X/F \right|, \, \zeta = ', num2str(z(3)),' $'],'Stiffness Line','Mass Line'},'Location','southwest','fontsize',lgndsize)
grid on; grid minor;
axis([1 1e2 1e-6 1e-1])
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
subplot(1,2,2)
plot(f,rad2deg(angle(Hjw(1,:))),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,rad2deg(angle(Hjw(2,:))),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,rad2deg(angle(Hjw(3,:))),'b','linewidth', 2), hold on % Zeta = 0.001
plot(f,rad2deg(-pi/2)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
xlabel('$f$ [Hz]')
ylabel('Receptance $\phi$ [$^{\circ}$]')
legend(['$\angle \left( X/F \right), \, \zeta = ', num2str(z(1)),' $'],['$\angle \left( X/F \right), \, \zeta = ', num2str(z(2)),' $'],['$\angle \left( X/F \right), \, \zeta = ', num2str(z(3)),' $'],'Location','northeast','fontsize',lgndsize)
grid on, grid minor, xlim([0 50]), ylim([-185 5])
set(gca,'fontsize',txtsize,'Ytick',[-180 -135 -90 -45 0],'XColor','k','YColor','k','ZColor','k','GridColor','k')

% REAL(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,real(Hjw(st,:)),'m','linewidth', 2)
        axis([0 30 -0.3e-3 +0.3e-3])
    elseif st == 2
        plot(f,real(Hjw(st,:)),'k','linewidth', 2)
        axis([0 30 -0.3e-2 +0.3e-2])
    else
        plot(f,real(Hjw(st,:)),'b','linewidth', 2)
        axis([0 30 -0.3e-1 +0.3e-1])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{X}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Re\{\frac{X}{F}\}$','Location','southwest','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% IM(X/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,imag(Hjw(st,:)),'m','linewidth', 2)
        axis([0 30 -0.00055 0])
        set(gca,'YTick',[-0.00055:0.0001:0]) %#ok<*NBRAK>
    elseif st == 2
        plot(f,imag(Hjw(st,:)),'k','linewidth', 2)
        axis([0 30 -0.0055 0])
        set(gca,'YTick',[-0.0055:0.001:0])
    else
        plot(f,imag(Hjw(st,:)),'b','linewidth', 2)
        axis([0 30 -0.055 0])
        set(gca,'YTick',[-0.055:0.01:0])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Im\{\frac{X}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Im\{\frac{X}{F}\}$','Location','southwest','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% REAL(X/F) PELA FREQUENCIA E IM(X/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
plot(f,real(Hjw(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,real(Hjw(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,real(Hjw(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
xlabel('$f$ [Hz]')
ylabel('$Re\{\frac{X}{F}\}$')
legend(['$Re\{X/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Re\{X/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Re\{X/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','southeast','fontsize',lgndsize)
grid on, grid minor,
axis([14 18 -0.3e-1 +0.3e-1])
set(gca,'fontsize',txtsize)
subplot(1,2,2)
plot(f,imag(Hjw(1,:)),'m','linewidth', 2), hold on
plot(f,imag(Hjw(2,:)),'k','linewidth', 2), hold on
plot(f,imag(Hjw(3,:)),'b','linewidth', 2), hold on
axis([14 18 -0.051 0])
xlabel('$f$ [Hz]')
ylabel('$Im\{\frac{X}{F}\}$')
legend(['$Im\{X/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Im\{X/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Im\{X/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','southeast','fontsize',lgndsize)
grid on, grid minor
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')

% NYQUIST - SEPARADOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(real(Hjw(st,:)),imag(Hjw(st,:)),'m','linewidth', 2)
    elseif st==2
        plot(real(Hjw(st,:)),imag(Hjw(st,:)),'k','linewidth', 2)
    else
        plot(real(Hjw(st,:)),imag(Hjw(st,:)),'b','linewidth', 2)
    end
    xlabel('$Re\{\frac{X}{F}\}$')
    ylabel('$Im\{\frac{X}{F}\}$')
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal') 
    grid on, grid minor
    axis image
end

% NYQUIST - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
plot(real(Hjw(1,:)),imag(Hjw(1,:)),'m','linewidth', 2),hold on
plot(real(Hjw(2,:)),imag(Hjw(2,:)),'k','linewidth', 2),hold on
plot(real(Hjw(3,:)),imag(Hjw(3,:)),'b','linewidth', 2),hold on
xlabel('$Re\{\frac{X}{F}\}$')
ylabel('$Im\{\frac{X}{F}\}$')
legend(['$\zeta = ', num2str(z(1)),' $'],['$\zeta = ', num2str(z(2)),' $'],['$\zeta = ', num2str(z(3)),' $'],'Location','southeast','fontsize',lgndsize)
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
grid on, grid minor
axis image

% NYQUIST - 3D
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot3(f,real(Hjw(1,:)),imag(Hjw(1,:)),'m','linewidth', 3),hold on
        plot3(f,real(Hjw(1,:)),-0.00055.*ones(1,length(w)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        plot3(f,-0.0003.*ones(1,length(w)),imag(Hjw(1,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Hjw(1,:)),imag(Hjw(1,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -0.0003 0.0003 -0.00055 0])
    elseif st==2
        plot3(f,real(Hjw(2,:)),imag(Hjw(2,:)),'k','linewidth', 3),hold on
        plot3(f,real(Hjw(2,:)),-0.0055.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,-0.003.*ones(1,length(w)),imag(Hjw(2,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Hjw(2,:)),imag(Hjw(2,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -0.003 0.003 -0.0055 0])
    else
        plot3(f,real(Hjw(3,:)),imag(Hjw(3,:)),'b','linewidth', 3),hold on
        plot3(f,real(Hjw(3,:)),-0.055.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,-0.03.*ones(1,length(w)),imag(Hjw(3,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Hjw(3,:)),imag(Hjw(3,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -0.03 0.03 -0.055 0])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{X}{F}\}$')
    zlabel('$Im\{\frac{X}{F}\}$')
    set(gca,'fontsize',txtsize-5,'YDir','reverse','XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    grid on, grid minor
end

%% MODILIDADE
Hm = 1i*w.*Hjw;
% ABS(\dot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Hm(st,:)),'m','linewidth', 2), hold on
        axis([1e-1 1e3 1e-5 1e-1])
        set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3],'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
    elseif st==2
        loglog(f,abs(Hm(st,:)),'k','linewidth', 2), hold on
        axis([1e-1 1e3 1e-5 1e0])
        set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3],'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
    else
        loglog(f,abs(Hm(st,:)),'b','linewidth', 2), hold on
        axis([1e-1 1e3 1e-5 1e1])
        set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3],'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])
    end
    loglog(f,abs(1i*w)/k.*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Stiffness
    loglog(f,abs(1i*w)./(w.^2*m),'--','linewidth', 1.8), hold on % Mass
    xlabel('$f$ [Hz]')
    ylabel('Mobility [m/N.s]')
    legend('$\left| \dot{X}/F \right|$','Stiffness Line','Mass Line','Location','southeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ANGULO(\dot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(f,rad2deg(angle(Hm(st,:))),'m','linewidth', 2), hold on
    elseif st==2
        plot(f,rad2deg(angle(Hm(st,:))),'k','linewidth', 2), hold on
    else
        plot(f,rad2deg(angle(Hm(st,:))),'b','linewidth', 2), hold on
    end
    plot(f,rad2deg(0)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
    xlabel('$f$ [Hz]')
    ylabel('Mobility $\phi$ [$^{\circ}$]')
    legend('$\angle{\left(\dot{X}/F\right)}$','Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), 
    grid on, grid minor, xlim([0 50]), ylim([-95 95])
    set(gca,'fontsize',txtsize,'Ytick',[-90 -45 0 45 90],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ABS(\dot{X}/F) PELA FREQUENCIA E ANGULO(\dot{X}/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
loglog(f,abs(Hm(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
loglog(f,abs(Hm(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
loglog(f,abs(Hm(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
loglog(f,abs(1i*w)/k.*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Rigidez
loglog(f,abs(1i*w)./(w.^2*m),'--','linewidth', 1.8), hold on % Massa
xlabel('$f$ [Hz]')
ylabel('Mobility [m/N.s]')
legend({['$\left| \dot{X}/F \right|, \, \zeta = ', num2str(z(1)),' $'],['$\left| \dot{X}/F \right|, \, \zeta = ', num2str(z(2)),' $'],['$\left| \dot{X}/F \right|, \, \zeta = ', num2str(z(3)),' $'],'Stiffness Line','Mass Line'},'Location','northeast','fontsize',lgndsize)
grid on; grid minor;
axis([1e-1 1e3 1e-4 1e1])
set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3],'YTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1])
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
subplot(1,2,2)
plot(f,rad2deg(angle(Hm(1,:))),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,rad2deg(angle(Hm(2,:))),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,rad2deg(angle(Hm(3,:))),'b','linewidth', 2), hold on % Zeta = 0.001
plot(f,rad2deg(0)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
xlabel('$f$ [Hz]')
ylabel('Mobility $\phi$ [$^{\circ}$]')
legend(['$\angle{\left( \dot{X}/F \right)}, \, \zeta = ', num2str(z(1)),' $'],['$\angle \left( \dot{X}/F \right), \, \zeta = ', num2str(z(2)),' $'],['$\angle \left( \dot{X}/F \right), \, \zeta = ', num2str(z(3)),' $'],'Location','northeast','fontsize',lgndsize)
grid on, grid minor, xlim([0 50]), ylim([-95 95])
set(gca,'fontsize',txtsize,'Ytick',[-90 -45 0 45 90],'XColor','k','YColor','k','ZColor','k','GridColor','k')

% REAL(\dot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,real(Hm(st,:)),'m','linewidth', 2)
        axis([0 40 0 0.055])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[0:0.01:0.055])
    elseif st == 2
        plot(f,real(Hm(st,:)),'k','linewidth', 2)
        axis([0 40 0 0.55])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[0:0.1:0.55])
    else
        plot(f,real(Hm(st,:)),'b','linewidth', 2)
        axis([0 40 0 5.5])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[0:1:5.5])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{\dot{X}}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Re\{\frac{\dot{X}}{F}\}$','Location','northeast','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% IM(\dot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,imag(Hm(st,:)),'m','linewidth', 2)
        axis([0 40 -0.03 0.03])
        set(gca,'YTick',[-0.03:0.01:0.03])
    elseif st == 2
        plot(f,imag(Hm(st,:)),'k','linewidth', 2)
        axis([0 40 -0.3 0.3])
        set(gca,'YTick',[-0.3:0.1:0.3])
    else
        plot(f,imag(Hm(st,:)),'b','linewidth', 2)
        axis([0 40 -3 3])
        set(gca,'YTick',[-3:1:3])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Im\{\frac{\dot{X}}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Im\{\frac{\dot{X}}{F}\}$','Location','northeast','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% REAL(\dot{X}/F) PELA FREQUENCIA E IM(\dot{X}/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
plot(f,real(Hm(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,real(Hm(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,real(Hm(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
xlabel('$f$ [Hz]')
ylabel('$Re\{\frac{\dot{X}}{F}\}$')
legend(['$Re\{\dot{X}/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Re\{\dot{X}/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Re\{\dot{X}/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','northeast','fontsize',lgndsize)
grid on, grid minor,
set(gca,'YTick',[0:0.5:5.5],'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
axis([12 20 0 5.5])
subplot(1,2,2)
plot(f,imag(Hm(1,:)),'m','linewidth', 2), hold on
plot(f,imag(Hm(2,:)),'k','linewidth', 2), hold on
plot(f,imag(Hm(3,:)),'b','linewidth', 2), hold on
axis([12 20 -3 3])
set(gca,'YTick',[-3:0.5:3],'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
xlabel('$f$ [Hz]')
ylabel('$Im\{\frac{\dot{X}}{F}\}$')
legend(['$Im\{\dot{X}/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Im\{\dot{X}/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Im\{\dot{X}/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','northeast','fontsize',lgndsize)
grid on, grid minor

% NYQUIST - SEPARADOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(real(Hm(st,:)),imag(Hm(st,:)),'m','linewidth', 2)
    elseif st==2
        plot(real(Hm(st,:)),imag(Hm(st,:)),'k','linewidth', 2)
    else
        plot(real(Hm(st,:)),imag(Hm(st,:)),'b','linewidth', 2)
    end
    xlabel('$Re\{\frac{\dot{X}}{F}\}$')
    ylabel('$Im\{\frac{\dot{X}}{F}\}$')
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal') 
    grid on, grid minor
    axis image
end

% NYQUIST - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
plot(real(Hm(1,:)),imag(Hm(1,:)),'m','linewidth', 2),hold on
plot(real(Hm(2,:)),imag(Hm(2,:)),'k','linewidth', 2),hold on
plot(real(Hm(3,:)),imag(Hm(3,:)),'b','linewidth', 2),hold on
xlabel('$Re\{\frac{\dot{X}}{F}\}$')
ylabel('$Im\{\frac{\dot{X}}{F}\}$')
legend(['$\zeta = ', num2str(z(1)),' $'],['$\zeta = ', num2str(z(2)),' $'],['$\zeta = ', num2str(z(3)),' $'],'Location','southeast','fontsize',lgndsize)
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
grid on, grid minor
axis image

% NYQUIST - 3D
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot3(f,real(Hm(1,:)),imag(Hm(1,:)),'m','linewidth', 3),hold on
        plot3(f,real(Hm(1,:)),-0.03.*ones(1,length(w)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        plot3(f,0.*ones(1,length(w)),imag(Hm(1,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        plot3(40.*ones(1,length(w)),real(Hm(1,:)),imag(Hm(1,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 40 0 0.05 -0.03 0.03])
    elseif st==2
        plot3(f,real(Hm(2,:)),imag(Hm(2,:)),'k','linewidth', 3),hold on
        plot3(f,real(Hm(2,:)),-0.3.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,0.*ones(1,length(w)),imag(Hm(2,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(40.*ones(1,length(w)),real(Hm(2,:)),imag(Hm(2,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 40 0 0.5 -0.3 0.3])
    else
        plot3(f,real(Hm(3,:)),imag(Hm(3,:)),'b','linewidth', 3),hold on
        plot3(f,real(Hm(3,:)),-3.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,0.*ones(1,length(w)),imag(Hm(3,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(40.*ones(1,length(w)),real(Hm(3,:)),imag(Hm(3,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 40 0 5 -3 3])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{\dot{X}}{F}\}$')
    zlabel('$Im\{\frac{\dot{X}}{F}\}$')
    set(gca,'fontsize',txtsize-5,'YDir','reverse','XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    grid on, grid minor
end

%% ACELERANCIA
Ha = -w.^2.*Hjw;
% ABS(\ddot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        loglog(f,abs(Ha(st,:)),'m','linewidth', 2), hold on
        axis([1e0 1e3 1e-3 1e1])
        set(gca,'XTick',[1e0 1e1 1e2 1e3],'YTick',[1e-3 1e-2 1e-1 1e0 1e1])
    elseif st==2
        loglog(f,abs(Ha(st,:)),'k','linewidth', 2), hold on
        axis([1e0 1e3 1e-3 1e2])
        set(gca,'XTick',[1e0 1e1 1e2 1e3],'YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2])
    else
        loglog(f,abs(Ha(st,:)),'b','linewidth', 2), hold on
        axis([1e0 1e3 1e-3 1e3])
        set(gca,'XTick',[1e0 1e1 1e2 1e3],'YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
    end
    loglog(f,abs((1i*w).^2)/k.*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Rigidez
    loglog(f,abs((1i*w).^2)./(w.^2*m),'--','linewidth', 1.8), hold on % Massa
    xlabel('$f$ [Hz]')
    ylabel('Accelerance [m/N.s$^2$]')
    legend('$\left| \ddot{X}/F \right|$','Stiffness Line','Mass Line','Location','southeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ANGULO(\ddot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(f,rad2deg(angle(Ha(st,:))),'m','linewidth', 2), hold on
    elseif st==2
        plot(f,rad2deg(angle(Ha(st,:))),'k','linewidth', 2), hold on
    else
        plot(f,rad2deg(angle(Ha(st,:))),'b','linewidth', 2), hold on
    end
    plot(f,rad2deg(pi/2)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
    xlabel('$f$ [Hz]')
    ylabel('Accelerance $\phi$ [$^{\circ}$]')
    legend('$\angle{\left(\ddot{X}/F\right)}$','Location','northeast','fontsize',lgndsize)
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal'), 
    grid on, grid minor, xlim([0 50]), ylim([-5 185])
    set(gca,'fontsize',txtsize,'Ytick',[0 45 90 135 180],'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% ABS(\ddot{X}/F) PELA FREQUENCIA E ANGULO(\ddot{X}/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
loglog(f,abs(Ha(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
loglog(f,abs(Ha(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
loglog(f,abs(Ha(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
loglog(f,abs((1i*w).^2)/k.*ones(1,length(w)),'--r','linewidth', 1.8), hold on % Rigidez
loglog(f,abs((1i*w).^2)./(w.^2*m),'--','linewidth', 1.8), hold on % Massa
xlabel('$f$ [Hz]')
ylabel('Accelerance [m/N.s$^2$]')
legend({['$\left| \ddot{X}/F \right|, \, \zeta = ', num2str(z(1)),' $'],['$\left| \ddot{X}/F \right|, \, \zeta = ', num2str(z(2)),' $'],['$\left| \ddot{X}/F \right|, \, \zeta = ', num2str(z(3)),' $'],'Stiffness Line','Mass Line'},'Location','southeast','fontsize',lgndsize)
grid on; grid minor;
axis([1e0 1e3 1e-3 1e3])
set(gca,'XTick',[1e0 1e1 1e2 1e3],'YTick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
subplot(1,2,2)
plot(f,rad2deg(angle(Ha(1,:))),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,rad2deg(angle(Ha(2,:))),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,rad2deg(angle(Ha(3,:))),'b','linewidth', 2), hold on % Zeta = 0.001
plot(f,rad2deg(pi/2)*ones(1,length(w)),'--r','linewidth', 1.8), hold on
xlabel('$f$ [Hz]')
ylabel('Accelerance $\phi$ [$^{\circ}$]')
legend(['$\angle{\left( \ddot{X}/F \right)}, \, \zeta = ', num2str(z(1)),' $'],['$\angle \left( \ddot{X}/F \right), \, \zeta = ', num2str(z(2)),' $'],['$\angle \left( \ddot{X}/F \right), \, \zeta = ', num2str(z(3)),' $'],'Location','northeast','fontsize',lgndsize)
grid on, grid minor, xlim([0 50]), ylim([-5 185])
set(gca,'fontsize',txtsize,'Ytick',[0 45 90 135 180],'XColor','k','YColor','k','ZColor','k','GridColor','k')

% REAL(\ddot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,real(Ha(st,:)),'m','linewidth', 2)
        axis([0 40 -3 3])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[-3:1:3])
    elseif st == 2
        plot(f,real(Ha(st,:)),'k','linewidth', 2)
        axis([0 40 -30 30])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[-30:10:30])
    else
        plot(f,real(Ha(st,:)),'b','linewidth', 2)
        axis([0 40 -300 300])
        set(gca,'XTick',[0 10 20 30 40],'YTick',[-300:100:300])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{\ddot{X}}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Re\{\frac{\ddot{X}}{F}\}$','Location','southeast','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% IM(\ddot{X}/F) PELA FREQUENCIA
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st == 1
        plot(f,imag(Ha(st,:)),'m','linewidth', 2)
        axis([0 40 -0.2 5.2])
        set(gca,'YTick',[0:1:5])
    elseif st == 2
        plot(f,imag(Ha(st,:)),'k','linewidth', 2)
        axis([0 40 -2 52])
        set(gca,'YTick',[0:10:50])
    else
        plot(f,imag(Ha(st,:)),'b','linewidth', 2)
        axis([0 40 -20 520])
        set(gca,'YTick',[0:100:500])
    end
    xlabel('$f$ [Hz]','Interpreter','Latex')
    ylabel('$Im\{\frac{\ddot{X}}{F}\}$')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    legend('$Im\{\frac{\ddot{X}}{F}\}$','Location','northeast','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end

% REAL(\ddot{X}/F) PELA FREQUENCIA E IM(\ddot{X}/F) PELA FREQUENCIA - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
subplot(1,2,1)
plot(f,real(Ha(1,:)),'m','linewidth', 2), hold on % Zeta = 0.1
plot(f,real(Ha(2,:)),'k','linewidth', 2), hold on % Zeta = 0.01
plot(f,real(Ha(3,:)),'b','linewidth', 2), hold on % Zeta = 0.001
xlabel('$f$ [Hz]')
ylabel('$Re\{\frac{\ddot{X}}{F}\}$')
legend(['$Re\{\ddot{X}/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Re\{\ddot{X}/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Re\{\ddot{X}/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','northwest','fontsize',lgndsize)
grid on, grid minor,
set(gca,'YTick',[-300:50:300],'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
axis([0 30 -300 300])
subplot(1,2,2)
plot(f,imag(Ha(1,:)),'m','linewidth', 2), hold on
plot(f,imag(Ha(2,:)),'k','linewidth', 2), hold on
plot(f,imag(Ha(3,:)),'b','linewidth', 2), hold on
axis([0 30 -20 520])
set(gca,'YTick',[0:50:500],'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
xlabel('$f$ [Hz]')
ylabel('$Im\{\frac{\ddot{X}}{F}\}$')
legend(['$Im\{\ddot{X}/F\}, \, \zeta = ', num2str(z(1)),' $'],['$Im\{\ddot{X}/F\}, \, \zeta = ', num2str(z(2)),' $'],['$Im\{\ddot{X}/F\}, \, \zeta = ', num2str(z(3)),' $'],'Location','northwest','fontsize',lgndsize)
grid on, grid minor

% NYQUIST - SEPARADOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot(real(Ha(st,:)),imag(Ha(st,:)),'m','linewidth', 2)
    elseif st==2
        plot(real(Ha(st,:)),imag(Ha(st,:)),'k','linewidth', 2)
    else
        plot(real(Ha(st,:)),imag(Ha(st,:)),'b','linewidth', 2)
    end
    xlabel('$Re\{\frac{\ddot{X}}{F}\}$')
    ylabel('$Im\{\frac{\ddot{X}}{F}\}$')
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal') 
    grid on, grid minor
    axis image
end

% NYQUIST - JUNTOS
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1])
plot(real(Ha(1,:)),imag(Ha(1,:)),'m','linewidth', 2),hold on
plot(real(Ha(2,:)),imag(Ha(2,:)),'k','linewidth', 2),hold on
plot(real(Ha(3,:)),imag(Ha(3,:)),'b','linewidth', 2),hold on
xlabel('$Re\{\frac{\ddot{X}}{F}\}$')
ylabel('$Im\{\frac{\ddot{X}}{F}\}$')
legend(['$\zeta = ', num2str(z(1)),' $'],['$\zeta = ', num2str(z(2)),' $'],['$\zeta = ', num2str(z(3)),' $'],'Location','southeast','fontsize',lgndsize)
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
grid on, grid minor
axis image

% NYQUIST - 3D
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.7])
for st=1:3
    subplot(1,3,st)
    if st==1
        plot3(f,real(Ha(1,:)),imag(Ha(1,:)),'m','linewidth', 3),hold on
        plot3(f,real(Ha(1,:)),0.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,-3.*ones(1,length(w)),imag(Ha(1,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Ha(1,:)),imag(Ha(1,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -3 3 0 5.1])
    elseif st==2
        plot3(f,real(Ha(2,:)),imag(Ha(2,:)),'k','linewidth', 3),hold on
        plot3(f,real(Ha(2,:)),0.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,-30.*ones(1,length(w)),imag(Ha(2,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Ha(2,:)),imag(Ha(2,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -30 30 0 51])
    else
        plot3(f,real(Ha(3,:)),imag(Ha(3,:)),'b','linewidth', 3),hold on
        plot3(f,real(Ha(3,:)),0.*ones(1,length(w)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(f,-300.*ones(1,length(w)),imag(Ha(3,:)),'-','linewidth', 1,'Color',[0 0 0]+0.7),hold on
        plot3(30.*ones(1,length(w)),real(Ha(3,:)),imag(Ha(3,:)),'-','linewidth', 0.5,'Color',[0 0 0]+0.7),hold on
        axis([0 30 -300 300 0 510])
    end
    xlabel('$f$ [Hz]')
    ylabel('$Re\{\frac{\ddot{X}}{F}\}$')
    zlabel('$Im\{\frac{\ddot{X}}{F}\}$')
    set(gca,'fontsize',txtsize-5,'YDir','reverse','XColor','k','YColor','k','ZColor','k','GridColor','k')
    title(['$\zeta = ', num2str(z(st)),' $'],'FontWeight','normal')
    grid on, grid minor
end