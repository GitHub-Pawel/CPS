clear all;
close all;
clc;

%% Dane
f_s = 256000;   %czêstotliwoœæ próbkowania w Hz
N = 2560;       

%% Filtr LP Butterworth
% Wyznaczam
[b,a] = butter(7,0.24);   %(kolejnosc filtrowania, czestotliwosc odciecia)
H=freqz(b, a, f_s);

% Rysuje wykres
figure(1);
subplot(2,1,1);
hold all;
plot(mag2db(abs(H)),'k');
plot([f_s/2 f_s/2], [min(mag2db(abs(H))) -40], 'g');
plot([0 f_s/4], [-3 -3], 'r');
title('Filtr Butterwortha');
legend('butter(...)');
xlabel('\omega [Hz]');
ylabel('[dB]')
subplot(2,1,2);
plot(b,'r o');
title('Rozklad biegunow');
xlabel('Real part');
ylabel('Imaginary part');

%% Filtr LP Czebyszew 1
% Wyznaczam
[b,a] = cheby1(4, 3, 0.29);   %(kolejnosc filtrowania,P2P têtnienie pasma przepustowego, czestotliwosc graniczna pasma)
H=freqz(b,a,f_s);

% Rysuje wykres
figure(2);
subplot(2,1,1);
hold all;
plot(mag2db(abs(H)),'k');
plot([f_s/2 f_s/2], [min(mag2db(abs(H))) -40], 'g');
plot([0 f_s/4], [-3 -3], 'r');
title('Filtr Czebyszewa 1');
legend('cheby1(...)');
xlabel('\omega [Hz]');
ylabel('[dB]')
subplot(2,1,2);
plot(b,'r o');
title('Rozklad biegunow');
xlabel('Real part');
ylabel('Imaginary part');

%%  Filtr LP Czebyszew 2
% Wyznaczam
[b,a] = cheby2(6, 40, 0.5);    %(kolejnosc filtrowania,P2P têtnienie pasma przepustowego, czestotliwosc graniczna pasma)
H=freqz(b, a, f_s);

% Rysuje wykres
figure(3);
subplot(2,1,1);
hold all;
plot(mag2db(abs(H)),'k');
plot([f_s/2 f_s/2], [min(mag2db(abs(H))) -40], 'g');
plot([0 f_s/4], [-3 -3], 'r');
title('Filtr Czebyszewa 2');
legend('cheby2(...)');
xlabel('\omega [Hz]');
ylabel('[dB]')
subplot(2,1,2);
plot(b,'r o');
title('Rozklad biegunow');
xlabel('Real part');
ylabel('Imaginary part');

%% Filtr LP Eliptyczny
% Wyznaczam
[b,a] = ellip(4, 3, 40, 0.4075);
H=freqz(b,a,f_s);

% Rysuje wykres
figure(4);
subplot(2,1,1);
hold all;
plot(mag2db(abs(H)),'k');
plot([f_s/2 f_s/2], [min(mag2db(abs(H))) -40], 'g');
plot([0 f_s/4], [-3 -3], 'r');
title('Filtr Elliptyczny');
legend('ellip(...)');
xlabel('\omega [Hz]');
ylabel('[dB]')
subplot(2,1,2);
plot(b,'r o');
title('Rozklad biegunow');
xlabel('Real part');
ylabel('Imaginary part');
