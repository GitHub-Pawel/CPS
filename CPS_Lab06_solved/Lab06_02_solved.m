clear all;
close all; 
clc;

[s9, fs] = audioread('s9.wav');

%% Rysuje spektogram
figure(1);
spectrogram(s9, 4096, 4096-512, [0:5:2000], fs);
hold all;
plot([0.697 0.697], [0 6], 'r--');
plot([0.770 0.770], [0 6], 'b--');
plot([0.852 0.852], [0 6], 'g--');
plot([0.941 0.941], [0 6], 'c--');
plot([1.209 1.209], [0 6], 'm-.');
plot([1.336 1.336], [0 6], 'w-.');
plot([1.447 1.447], [0 6], 'k-.');
legend('697 Hz', '770 Hz', '852 Hz', '941 Hz', '1209 Hz', '1336 Hz', '1447 Hz');

%% Rozkodowany PIN
%   9   1   5   2   0






%% Kod z zadania 1
f_s = 16000;                %[Hz]
n = f_s;                    %liczba prboek
f_down = 1189;              %dolna czestotliwosc graniczna [Hz]
f_up =  1229;               %gorna czestotliwosc graniczna [Hz]
f = 0:n-1;       
w = f*2*pi;               %dziedzina [s]
B_w = f_up-f_down;          %szerokosc pasma
f_c = (f_up + f_down)/2;    %czestotliwosc srodkowa
t = 0:1/f_s:(length(s9)-1)/f_s;  %dziedzina sygnalu x [s]

load('butter.mat');
% z - zera
% p - bieguny
% k - wspolczynnik wzmocnienia analogowego filtru

%********************************
[b_s, a_s] = zp2tf(z, p, k); 
H_s = freqs(b_s, a_s, w);  
%********************************


%% Transformata biliniowa
%*************************************
[b_d, a_d] = bilinear(b_s, a_s, f_s);
H_d = freqz(b_d, a_d, f, f_s); 
%*************************************

%% Filtrowanie sygnalu
s9_filtered = filter(b_d, a_d, s9)*2.5144;
figure(2);
spectrogram(s9_filtered, 4096, 4096-512, [0:5:2000], fs);
hold all;
plot([0.697 0.697], [0 6], 'r--');
plot([0.770 0.770], [0 6], 'b--');
plot([0.852 0.852], [0 6], 'g--');
plot([0.941 0.941], [0 6], 'c--');
plot([1.209 1.209], [0 6], 'm-.');
plot([1.336 1.336], [0 6], 'w-.');
plot([1.447 1.447], [0 6], 'k-.');
legend('697 Hz', '770 Hz', '852 Hz', '941 Hz', '1209 Hz', '1336 Hz', '1447 Hz');


%% Wykresy sygnalow

figure(3);
hold all;
plot(t, s9, 'r--');
plot(t(1:(length(s9_filtered))-311), s9_filtered(312:length(s9_filtered)), 'c-.');
title('Porownanie sygnalow przed i po filtracji');
legend('s9 przed filtracja', 's9 po filtracji');
xlabel('czas [s]');
grid;