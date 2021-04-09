%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK #2
% Docente: Michael John Brennan
% Discente: Estevão Fuzaro de Almeida
% Data: 18/03/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% VARIAVEIS
m = 1; %[kg]
k = 1e4; %[N/m]
z = [0.001 0.01 0.1]; %[adimensional]
t = 0:0.001:100; %[s]
f = 0:0.0001:60; %[Hz]
w = 2*pi*f; %[rad/s]
figc = 1; % contador de figuras

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m); %[rad/s]
fn = wn/(2*pi); %[Hz]
wd = wn*sqrt(1-z.^2); %[rad/s]
c = z*2*sqrt(k*m); %[Ns/m]

%% 1º METODO PARA ESTIMAR AMORTECIMENTO - DOMINIO DO TEMPO
% OBTENCAO DE h(t) E DO ENVELOPE
ht = []; env = []; % Criando os vetores
for st=1:3
    ht = [ht; 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t)]; %#ok<*AGROW>
    env = [env; 1/(m*wd(st))*exp(-z(st)*wn*t)];
end
% APLICANDO TRANSFORMADA DE HILBERT
figure(figc); figc = figc + 1;
for st=1:3
    subplot(2,3,st)
    X = hilbert(ht(st,:)); % HILBERT
    if st==1
        tp = t(5000:45000); % PARTE DO TEMPO
        Xp = abs(X(5000:45000)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'m','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])
    elseif st==2
        tp = t(500:4500); % PARTE DO TEMPO
        Xp = abs(X(500:4500)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'k','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])
    else
        tp = t(50:320); % PARTE DO TEMPO
        Xp = abs(X(50:320)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'b','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])    
    end
    xlabel('$t$ [s]')
    ylabel('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal');
    legend('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$','Part of data')
    grid on, grid minor
    set(gca,'fontsize',18)
end
for st=1:3
    subplot(2,3,st+3)
    X = hilbert(ht(st,:)); % HILBERT
    if st==1
        tp = t(5000:45000); % PARTE DO TEMPO
        Xp = abs(X(5000:45000)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'m','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 50])
    elseif st==2
        tp = t(500:4500); % PARTE DO TEMPO
        Xp = abs(X(500:4500)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'k','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 5])
    else
        tp = t(50:320); % PARTE DO TEMPO
        Xp = abs(X(50:320)); % PARTE DA TRANSFORMADA
        plot(t,log(abs(X)),'b','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 0.5])    
    end
    xlabel('$t$ [s]')
    ylabel('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal');
    legend('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$','Part of data')
    grid on, grid minor
    set(gca,'fontsize',18)
    % ESTIMANDO ZETA PELO TEMPO
        grad = polyfit(tp,log(Xp),1);
        z_est_time(st) = -grad(1)/wn; %#ok<*SAGROW>
end

%% 2º METODO PARA ESTIMAR AMORTECIMENTO - DOMINIO DA FREQUENCIA
% OBTENCAO DE H(jw)
Hjw = [];
for st=1:3
    Hjw = [Hjw; 1./(k-w.^2*m+1i*w*c(st))];
end
% ANGULO(X/F) PELA FREQUENCIA COM ZOOM
for st=1:3
    fig = figure(figc); figc = figc + 1;
    plot(f,rad2deg(-pi/4)*ones(1,length(w)),'r--','linewidth', 1.5), hold on
    plot(f,rad2deg(-pi/2)*ones(1,length(w)),'r--','linewidth', 1.5), hold on
    plot(f,rad2deg(-3*pi/4)*ones(1,length(w)),'r--','linewidth', 1.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('$\phi$ [$^{\circ}$]')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal'), 
    grid on, grid minor
    xlim([0 60])
    set(gca,'fontsize',18,'YTick',[-180 -135 -90 -45 0])
    if st==1
        plot(f,rad2deg(angle(Hjw(st,:))),'m','linewidth', 2.5), hold on
        xlm = [15.89,15.94];
        ylm = [-140,-40];
        pos = [0.50 0.20 0.38 0.635];
        annotation('arrow',[0.3555,0.5],[0.45,0.48]);
        magnifyPlot(xlm,ylm,pos,1.8,true,xlm(1):0.01:xlm(2),[-135 -120 -105 -90 -75 -60 -45])
            % ESTIMANDO ZETA PELA FREQUENCIA
            z_est_freq(st) = (15.9313-15.8995)/(2*fn);
    elseif st==2
        plot(f,rad2deg(angle(Hjw(st,:))),'k','linewidth', 2.5), hold on
        xlm = [15.7,16.1];
        ylm = [-140,-40];
        pos = [0.50 0.20 0.38 0.635];
        annotation('arrow',[0.3555,0.5],[0.45,0.48]);
        magnifyPlot(xlm,ylm,pos,1.8,true,xlm(1):0.1:xlm(2),[-135 -120 -105 -90 -75 -60 -45])
            % ESTIMANDO ZETA PELA FREQUENCIA
            z_est_freq(st) = (16.0748-15.7559)/(2*fn);
    else
        plot(f,rad2deg(angle(Hjw(st,:))),'b','linewidth', 2.5), hold on
        xlm = [14.3,17.9];
        ylm = [-140,-40];
        pos = [0.50 0.20 0.38 0.635];
        annotation('arrow',[0.3555,0.5],[0.45,0.48]);
        magnifyPlot(xlm,ylm,pos,1.8,true,xlm(1):0.4:xlm(2),[-135 -120 -105 -90 -75 -60 -45])
            % ESTIMANDO ZETA PELA FREQUENCIA
            z_est_freq(st) = (17.5755-14.3977)/(2*fn);
    end
end

%% COMPARACAO ENTRE OS METODOS DO TEMPO E DA FREQUENCIA
perc_error_time = abs(z_est_time - z)./z*100;
perc_error_freq = abs(z_est_freq - z)./z*100;

fprintf('Zeta estimados pelo tempo: %f \n', z_est_time)
fprintf('\n')
fprintf('Erros percentuais pelo tempo: %f \n', perc_error_time)
fprintf('\n')
fprintf('Zeta estimados pela frequencia: %f \n', z_est_freq)
fprintf('\n')
fprintf('Erros percentuais pela freq: %f \n', perc_error_freq)