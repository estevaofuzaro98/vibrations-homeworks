%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 2
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 18/03/2021

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
z = [0.001 0.01 0.1];   % Zeta [adimensional]
t = 0:0.001:100;        % Vetor de Tempo [s]
f = 0:0.0001:60;        % Frequencia [Hz]
w = 2*pi*f;             % Freq. Angular [rad/s]

%% PARAMETROS DO SISTEMA
wn = sqrt(k/m);         % Freq. Natural [rad/s]
fn = wn/(2*pi);         % Freq. Natural [Hz]
wd = wn*sqrt(1-z.^2);   % Freq. Nat. Amortecida [rad/s]
c = 2*z*sqrt(k*m);      % Amortecimento [N.s/m]

%% 1st: ESTIMANDO AMORTECIMENTO PELO DOMINIO DO TEMPO
% OBTENCAO DE h(t)
ht = []; % Criando os vetores
for st=1:3
    ht(st,:) = 1/(m*wd(st))*exp(-z(st)*wn*t).*sin(wd(st)*t); %#ok<*SAGROW>
end
% APLICANDO TRANSFORMADA DE HILBERT
figure; X = [];
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    X(st,:) = hilbert(ht(st,:));    % HILBERT
    if st==1
        tp = t(5000:45000);         % PARTE DO TEMPO
        Xp = X(st,5000:45000);      % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'m','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])
    elseif st==2
        tp = t(500:4500);           % PARTE DO TEMPO
        Xp = abs(X(st,500:4500));   % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'k','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])
    else
        tp = t(50:320);             % PARTE DO TEMPO
        Xp = X(st,50:320);          % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'b','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 100])
    end
    xlabel('$t$ [s]')
    ylabel('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal');
    legend('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$','Part of data','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
end
figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
for st=1:3
    subplot(1,3,st)
    if st==1
        tp = t(5000:45000);         % PARTE DO TEMPO
        Xp = X(st,5000:45000);      % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'m','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 50])
    elseif st==2
        tp = t(500:4500);           % PARTE DO TEMPO
        Xp = X(st,500:4500);        % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'k','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 5])
    else
        tp = t(50:320);             % PARTE DO TEMPO
        Xp = X(st,50:320);          % PARTE DA TRANSFORMADA
        plot(t,log(abs(X(st,:))),'b','linewidth', 2), hold on
        plot(tp,log(abs(Xp)),'r','linewidth', 2.2), hold on
        xlim([0 0.5])    
    end
    xlabel('$t$ [s]')
    ylabel('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal');
    legend('$\ln \left| \mathcal{H}\left[h(t)\right] \right|$','Part of data','fontsize',lgndsize)
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    % ESTIMANDO ZETA PELO TEMPO
        grad = polyfit(tp,log(Xp),1);
        zTempo(st) = -real(grad(1))/wn; %#ok<*SAGROW>
end

%% 2nd: ESTIMANDO AMORTECIMENTO PELO DOMINIO DA FREQUENCIA
% OBTENCAO DE H(jw)
Hjw = [];
for st=1:3
    Hjw = [Hjw; 1./(k-w.^2*m+1i*w*c(st))];
end
% ABS(X/F) PELA FREQUENCIA COM ZOOM
for st=1:3
    maximo = max(abs(Hjw(st,:)));
    maximo_sqrt2 = (max(abs(Hjw(st,:)))/sqrt(2));
    figure
    set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.8])
    plot(f,maximo*ones(1,length(w)),'r--','linewidth', 1.5), hold on
    plot(f,maximo_sqrt2*ones(1,length(w)),'r--','linewidth', 1.5), hold on
    xlabel('$f$ [Hz]')
    ylabel('$\left| \frac{X}{F} \right|$')
    title(['$\zeta = ', num2str(z(st)), ' $'],'FontWeight','normal'), 
    grid on, grid minor
    set(gca,'fontsize',txtsize,'XColor','k','YColor','k','ZColor','k','GridColor','k')
    if st==1
        xlim([15 19])
        plot(f,abs(Hjw(st,:)),'m','linewidth', 2.5), hold on
        xlm = [15.89,15.94];
        ylm = [0.95*maximo_sqrt2,1.05*maximo];
        pos = [0.50 0.21 0.38 0.635];
        annotation('arrow',[0.313,0.5],[0.65,0.67]);
        magnifyPlot(xlm,ylm,pos,2.5,true,xlm(1):0.01:xlm(2),[maximo_sqrt2 maximo],lgndsize)
            % ESTIMANDO ZETA PELA FREQUENCIA
            zFreq(st) = (15.9313-15.8995)/(2*fn); % ALTERAR AQUI!!!
    elseif st==2
        xlim([14 24])
        plot(f,abs(Hjw(st,:)),'k','linewidth', 2.5), hold on
        xlm = [15.7,16.1];
        ylm = [0.95*maximo_sqrt2,1.05*maximo];
        pos = [0.50 0.21 0.38 0.635];
        annotation('arrow',[0.293,0.5],[0.65,0.67]);
        magnifyPlot(xlm,ylm,pos,2.5,true,xlm(1):0.1:xlm(2),[maximo_sqrt2 maximo],lgndsize)
            % ESTIMANDO ZETA PELA FREQUENCIA
            zFreq(st) = (16.0725-15.7535)/(2*fn); % ALTERAR AQUI!!!
    else
        xlim([5 60])
        plot(f,abs(Hjw(st,:)),'b','linewidth', 2.5), hold on
        xlm = [13.85,17.45];
        ylm = [0.95*maximo_sqrt2,1.05*maximo];
        pos = [0.50 0.21 0.38 0.635];
        annotation('arrow',[0.306,0.5],[0.65,0.67]);
        magnifyPlot(xlm,ylm,pos,2.5,true,xlm(1):0.6:xlm(2),[maximo_sqrt2 maximo],lgndsize)
            % ESTIMANDO ZETA PELA FREQUENCIA
            zFreq(st) = (17.2821-14.0573)/(2*fn); % ALTERAR AQUI!!!
    end
end

%% COMPARACAO ENTRE OS METODOS
for st=1:3
    erroTempo(st) = (abs(zTempo(st) - z(st))/z(st))*100;
    erroFreq(st) = (abs(zFreq(st) - z(st))/z(st))*100;
end
fprintf('Zetas estimados pelo tempo: %f \n', zTempo)
fprintf('\n')
fprintf('Zetas estimados pela frequencia: %f \n', zFreq)
fprintf('\n')
fprintf('Erros percentuais pelo tempo: %f \n', erroTempo)
fprintf('\n')
fprintf('Erros percentuais pela frequencia: %f \n', erroFreq)
fprintf('\n')