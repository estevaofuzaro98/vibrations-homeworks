%% MECHANICAL VIBRATIONS (2021/1) - HOMEWORK 5
% Docente: Michael John Brennan
% Discente: Estevao Fuzaro de Almeida
% Data: 15/04/2021

% INICIALIZACAO
clc; clear all; close all; format long; %#ok<*CLALL>
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
txtsize = 26;
lgndsize = 18;

%% VARIAVEIS
Fs = 5000;  % Freq. de Amostragem [Hz]
dt = 1/Fs;  % Incremento de Tempo [s]
T = 10;     % Periodo [s]
t = 0:dt:T; % Vetor de Tempo [s]

%% SINAL AVALIADO: SENO
Fn = 2;                         % Freq. de Excitacao [Hz]
ASeno = 1;                      % Amplitude do Sinal
wn = 2*pi*Fn;                   % Freq. de Excitacao [rad/s]
xSeno = ASeno*sin(wn*t);        % Sinal Seno
XSeno = fft(xSeno)*dt;          % Transformada de Fourier
PSDSeno = conj(XSeno).*XSeno/T; % Calculo de PSD
NFSeno = round(length(XSeno)/2);% Numero de Pontos de Freq.
PSDSeno = 2*PSDSeno(1:NFSeno);  % PSD One-Sided
fSeno = linspace(0,Fs/2,NFSeno);% Vetor de Frequencia [Hz]

figure  % Plotando no Tempo
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,xSeno,'b','linewidth',3), hold on
grid on, grid minor
xlabel('$t$ [s]')
ylabel('$x(t)$')
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

figure  % Plotando na Frequencia
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(fSeno,PSDSeno,'r','linewidth',3), hold on
grid on, grid minor
xlim([0 8])
xlabel('Frequency [Hz]')
ylabel('PSD [m$^2$/Hz]')
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

areaSeno = trapz(PSDSeno)/T;    % Area sob Curva PSD do Seno
areaSenoSQRT = sqrt(areaSeno);  % Raiz da Area
RMSSeno = rms(xSeno);           % RMS do Seno
fprintf('sqrt(A_PSD) Seno: %f \n', areaSenoSQRT)
fprintf('RMS Seno: %f \n', RMSSeno)
fprintf('\n')

%% SINAL AVALIADO: CHIRP
f1Chirp = 10;                           % Freq. Minima Chirp
f2Chirp = 100;                          % Freq. Maxima Chirp
xChirp = chirp(t,f1Chirp,T,f2Chirp);    % Sinal Chirp
XChirp = fft(xChirp)*dt;                % Transformada de Fourier
PSDChirp = conj(XChirp).*XChirp/T;      % Calculo de PSD
NFChirp = round(length(XChirp)/2);      % Numero de Pontos de Freq.
PSDChirp = 2*PSDChirp(1:NFChirp);       % PSD One-Sided
fChirp = linspace(0,Fs/2,NFChirp);      % Vetor de Frequencia [Hz]

figure  % Plotando no Tempo
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,xChirp,'b','linewidth',3), hold on
grid on, grid minor
xlabel('$t$ [s]')
ylabel('$x(t)$')
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

figure  % Plotando na Frequencia
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(fChirp,PSDChirp,'r','linewidth',3), hold on
grid on, grid minor
xlim([0 110])
xlabel('Frequency [Hz]')
ylabel('PSD [m$^2$/Hz]')
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

areaChirp = trapz(PSDChirp)/T;    % Area sob Curva PSD do Chirp
areaChirpSQRT = sqrt(areaChirp);  % Raiz da Area
RMSChirp = rms(xChirp);           % RMS do Chirp
fprintf('sqrt(A_PSD) Chirp: %f \n', areaChirpSQRT)
fprintf('RMS Chirp: %f \n', RMSChirp)
fprintf('\n')

%% SINAL AVALIADO: RANDOM
xRand = randn(1,length(t));         % Sinal Randomico
FreqCorte = 800;                    % Freq. de Corte Filtro [Hz]
OrdFiltro = 8;                      % Ordem do Filtro Low-Pass
a = FreqCorte/(Fs/2);               % Razão de Freq. Wn
[B,A] = butter(OrdFiltro,a,'low');  % Filtro Low-Pass
xRandFiltro = filter(B,A,xRand);    % Sinal Filtrado

% MÉTODO 1
N = length(xRandFiltro);                % Tamanho Sinal
xDFT = fft(xRandFiltro);                % Transformada de Fourier
xDFT = xDFT(1:round(N/2));              % Dividindo vetor em 50%
PSDRand = (1/(Fs*N))*abs(xDFT).^2;      % Calculo de PSD usando FFT
PSDRand(2:end-1) = 2*PSDRand(2:end-1);  % PSD One-Sided
fRand = 0:Fs/length(xRandFiltro):Fs/2;  % Vetor de Frequencia [Hz]

% MÉTODO 2
% XRand = fft(xRandFiltro)*dt;    % Transformada de Fourier
% PSDRand = conj(XRand).*XRand/T; % Calculo de PSD
% NFRand = round(length(XRand)/2);% Numero de Pontos de Freq.
% PSDRand = 2*PSDRand(1:NFRand);  % PSD One-Sided
% fRand = linspace(0,Fs/2,NFRand);% Vetor de Frequencia [Hz]

% MÉTODO 3
% NFFT = max(t)*Fs; NOverLap = round(NFFT/2);
% PSDRand = pwelch(xRandFiltro,hann(NFFT),NOverLap,NFFT,Fs,'onesided');  % Calculo de PSD
% fRand = linspace(0,Fs/2,length(PSDRand));                   % Vetor de Frequencia [Hz]

figure  % Plotando no Tempo
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
plot(t,xRand,'b',t,xRandFiltro,'r','linewidth',2)
grid on, grid minor
xlabel('$t$ [s]')
ylabel('$x(t)$')
legend({'Random Signal','Random Signal Filtered'},'Location','northeast','fontsize',lgndsize)
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

figure  % Plotando na Frequencia
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 0.6])
semilogy(fRand,PSDRand,'r','linewidth',3), hold on
grid on, grid minor
xlabel('Frequency [Hz]')
ylabel('PSD [m$^2$/Hz]')
set(gca,'fontsize',txtsize,'XColor','k','YColor','k','GridColor','k')

areaRand = trapz(PSDRand)/T;    % Area sob Curva PSD do Rand
areaRandSQRT = sqrt(areaRand);  % Raiz da Area
RMSRand = rms(xRandFiltro);     % RMS do Rand
fprintf('sqrt(A_PSD) Rand: %f \n', areaRandSQRT)
fprintf('RMS Rand: %f \n', RMSRand)
fprintf('\n')