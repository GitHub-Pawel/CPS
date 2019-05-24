close all;
clear all;
clc;

%% 1) Wybor filtra (BP) i okreslenie jego pulsacji granicznych (czyli transmitancji):
f_pr = 1200; %Czestotliwosc probkowania [Hz]
d_f = 200;   %Szerokosc pasma przepustowego [Hz]
f_c = 300;   %Czestotliwosc srodkowa [Hz]
N = 128;     %Liczba probek filtru
%N = 129;     %Liczba probek filtru
M = (N-1)/2;

f_1 = f_c - d_f/2;
f_2 = f_c + d_f/2;
F_1 = f_1/f_pr; 
F_2 = f_2/f_pr; 
O_1 = 2*pi*F_1;
O_2 = 2*pi*F_2;

n = -150:150;

f = (0:length(n)-1)*2*f_pr/(2*length(n));  %przeskalowanie dziedziny [Hz]

if mod(length(n), 2) == 0
    n_c = length(n)/2;
else
    n_c = (length(n)+1)/2;
end

if mod(N, 2) == 0
    Ni = N;
else
    Ni = N-1;
end

%% 2) Wyznaczenie dyskretnej odpowiedzi filtra h_PH(n):
h_HP(1:length(n)) = 2*F_2*sin(O_2*n)./(O_2*n) - 2*F_1*sin(O_1*n)./(O_1*n);
h_HP(n_c) = 2*(F_2 - F_1);

%% 2.5) Wyznaczenie dyskretnych okien czasowych w(n):
% Prostok¹tne:
w_pro(1:length(n)) = 0;
w_pro((n_c - Ni/2) : (n_c + Ni/2)) = 1;

% Hanninga:
w_han(1:length(n)) = 0.5 + 0.5*cos(2*pi*n ./ (2*M+1));
w_han(1:(n_c - Ni/2)) = 0;
w_han((n_c + Ni/2):length(n)) = 0;

% Hamminga:
w_ham(1:length(n)) = 0.54 + 0.46*cos(2*pi*n ./ (2*M+1));
w_ham(1:(n_c - Ni/2)) = 0;
w_ham((n_c + Ni/2):length(n)) = 0;

% Blackmana:
w_bla(1:length(n)) = 0.42 + 0.5*cos(2*pi*n ./ (2*M+1)) + 0.08*cos(4*pi*n ./ (2*M+1));
w_bla(1:(n_c - Ni/2)) = 0;
w_bla((n_c + Ni/2):length(n)) = 0;

% Blackmana-Harrisa:
w_bh(1:length(n)) = 0.35875 + 0.48829*cos(2*pi*n ./ (2*M+1)) + 0.14128*cos(4*pi*n ./ (2*M+1)) + 0.01168*cos(6*pi*n ./ (2*M+1));
w_bh(1:(n_c - Ni/2)) = 0;
w_bh((n_c + Ni/2):length(n)) = 0;

%% 3) wymno¿enia obliczonego h(n) z wybran¹ funkcj¹ okna czasowego w(n):

h_w_pro = h_HP .* w_pro; % Prostok¹tne
h_w_han = h_HP .* w_han; % Hanninga
h_w_ham = h_HP .* w_ham; % Hamminga
h_w_bla = h_HP .* w_bla; % Blackmana
h_w_bh = h_HP .* w_bh;   % Blackmana-Harrisa


%% 4) przesuniêcia hw(n) w prawo o M próbek i pobrania 2M+1 próbek: 
f = f - 2*length(n);

%% Charakterystyki amplitudowo-czestotliwosciowe:
figure(1);
hold all;
plot(f, mag2db(abs(fft(h_w_pro))), 'b-');
plot(f, mag2db(abs(fft(h_w_han))), 'r-');
plot(f, mag2db(abs(fft(h_w_ham))), 'g-');
plot(f, mag2db(abs(fft(h_w_bla))), 'c-');
plot(f, mag2db(abs(fft(h_w_bh))), 'm-');

plot([-f_1 -f_1], [-180 20], 'k:');
plot([-f_2 -f_2], [-180 20], 'k:');
plot([f_1 f_1], [-180 20], 'k:');
plot([f_2 f_2], [-180 20], 'k:');

title('Charakterystyki amplitudowo-czestotliwosciowe filtrow');
legend('Prostokat', 'Hanning', 'Hamming', 'Blackman', 'Blackmana-Harrisa');
xlabel('[Hz]');
ylabel('[dB]');

%% Charakterystyki fazowo-czestotliwosciowe:
figure(2);
hold all;
plot(f, unwrap(angle(fft(h_w_pro))), 'b-');
plot(f, unwrap(angle(fft(h_w_han))), 'r-');
plot(f, unwrap(angle(fft(h_w_ham))), 'g-');
plot(f, unwrap(angle(fft(h_w_bla))), 'c-');
plot(f, unwrap(angle(fft(h_w_bh))), 'm-');

plot([-f_1 -f_1], [-800 100], 'k:');
plot([-f_2 -f_2], [-800 100], 'k:');
plot([f_1 f_1], [-800 100], 'k:');
plot([f_2 f_2], [-800 100], 'k:');
title('Charakterystyki fazowo-czestotliwosciowe filtrow');
legend('Prostokat', 'Hanning', 'Hamming', 'Blackman', 'Blackmana-Harrisa');
xlabel('[Hz]');





%% ************** Sygnal, ktory zostanie poddany f iltracji ***************
f1 = 50; f2 = 300; f3 = 413; %czestotliwosci sygnalow skladowych
A1 = 2500; A2 = 500; A3 = 5000; %amplitudy sygnalow skladowych

sin1 = A1*sin(2*pi * f1 * [0:1/f_pr:(length(n)-1)/f_pr]);
sin2 = A2*sin(2*pi * f2 * [0:1/f_pr:(length(n)-1)/f_pr]);
sin3 = A3*sin(2*pi * f3 * [0:1/f_pr:(length(n)-1)/f_pr]);

x = sin1 + sin2 + sin3;

%% ***************** Filtracja sygnalu x roznymi filtrami *****************
x_pro = ifft(fft(x).*fft(h_w_pro));
x_han = ifft(fft(x).*fft(h_w_han));
x_ham = ifft(fft(x).*fft(h_w_ham));
x_bla = ifft(fft(x).*fft(h_w_bla));
x_bh = ifft(fft(x).*fft(h_w_bh));

%% *************** Porownanie sygnalu przed i po filtracja ****************
%figure(3);
%hold all;
%plot(n, x, 'k-.');
%plot(n, x_pro, 'b-');
%title('Sygnal przed i po filtracji');
%legend('Sygnal przed filtracja', 'Sygnal po filtracji [prostokat]');
%xlabel('czas [s]');


%figure(4);
%hold all;
%plot(n, x, 'k-.');
%plot(n, x_han, 'r-');
%title('Sygnal przed i po filtracji');
%legend('Sygnal przed filtracja', 'Sygnal po filtracji [hanning]');
%xlabel('czas [s]');

figure(5);
hold all;
plot(n, x, 'k-.');
plot(n, x_ham, 'g-');
title('Sygnal przed i po filtracji');
legend('Sygnal przed filtracja', 'Sygnal po filtracji [hamming]');
xlabel('czas [s]');

%figure(6);
%hold all;
%plot(n, x, 'k-.');
%plot(n, x_bla, 'c-');
%title('Sygnal przed i po filtracji');
%legend('Sygnal przed filtracja', 'Sygnal po filtracji [blackman]');
%xlabel('czas [s]');

%figure(7);
%hold all;
%plot(n, x, 'k-.');
%plot(n, x_bh, 'm-');
%title('Sygnal przed i po filtracji');
%legend('Sygnal przed filtracja', 'Sygnal po filtracji [blackmana-harris]');
%xlabel('czas [s]');

%% ************************ Cha-ki a-cz sygnalu x *************************
figure(8);
hold all;
plot (f, abs(fft(x)), 'k:'); 
plot(f, abs(fft(x_pro)), 'b-o');
plot(f, abs(fft(x_han)), 'r-x');
plot(f, abs(fft(x_ham)), 'g-*');
plot(f, abs(fft(x_bla)), 'c-s');
plot(f, abs(fft(x_bh)), 'm-d');
title('Charakterystyki amplitudowo-czestotliwosciowe sygnalu x przefiltrowanego przez rozne filtry');
legend('X Przed filtracja', 'Prostokat', 'Hanning', 'Hamming', 'Blackman', 'Blackmana-Harrisa');
xlabel('[Hz]');
%ylabel('[dB]');

%% ************************ Cha-ki f-cz sygnalu x *************************
figure(9);
hold all;
plot(f, unwrap(angle(fft(x_pro))), 'b-');
plot(f, unwrap(angle(fft(x_han))), 'r-');
plot(f, unwrap(angle(fft(x_ham))), 'g-');
plot(f, unwrap(angle(fft(x_bla))), 'c-');
plot(f, unwrap(angle(fft(x_bh))), 'm-');

title('Charakterystyki fazowo-czestotliwosciowe  sygnalu x przefiltrowanego przez rozne filtry');
legend('Prostokat', 'Hanning', 'Hamming', 'Blackman', 'Blackmana-Harrisa');
xlabel('[Hz]');

%% ********* Widmowa gestosc mocy sygnalu x przed i po filtracja **********
[pxx_x, f_x] = periodogram(x, [], length(n), f_pr);
pxx_pro = periodogram(x_pro, [], length(n), f_pr);
pxx_han = periodogram(x_han, [], length(n), f_pr);
pxx_ham = periodogram(x_ham, [], length(n), f_pr);
pxx_bla = periodogram(x_bla, [], length(n), f_pr);
pxx_bh = periodogram(x_bh, [], length(n), f_pr);


figure(10);
hold all;
plot(f_x, mag2db(pxx_pro), 'b-o');
plot(f_x, mag2db(pxx_pro), 'r-x');
plot(f_x, mag2db(pxx_pro), 'g-*');
plot(f_x, mag2db(pxx_pro), 'c-s');
plot(f_x, mag2db(pxx_pro), 'm--');
plot(f_x, mag2db(pxx_x), 'k:');
title('Widmowa gestosc energi');
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)')
grid;
legend('Sygnal po filtracji (prostokat)', 'Sygnal po filtracji (hanning)', 'Sygnal po filtracji (hamming)', 'Sygnal po filtracji (blackman)', 'Sygnal po filtracji (blackman-harris)', 'Sygnal przed filtracja');

%% ************************************************************************