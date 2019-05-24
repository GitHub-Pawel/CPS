%BP - Bandpass - pasmowo-przepustowy
clear all;
close all;
clc;

f_s = 16000;                %[Hz]
n = f_s;                    %liczba prboek
f_down = 1189;              %dolna czestotliwosc graniczna [Hz]
f_up =  1229;               %gorna czestotliwosc graniczna [Hz]
f = 0:n-1;       
w = f*2*pi;               %dziedzina [s]
B_w = f_up-f_down;          %szerokosc pasma
f_c = (f_up + f_down)/2;    %czestotliwosc srodkowa

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

%% Rysuje wykresy

figure(1);
hold all;
plot(f, mag2db(abs(H_s)), 'b--');
plot(f, mag2db(abs(H_d)), 'r-.');
title('Charakterystyka amplitudowo-czestotliwosciowa transmitancji');
legend('Transmitancja pseudo analogowa', 'Transmitancja cyfrowa');
xlabel('[Hz]');
ylabel('[dB]');
grid;


graniczna_analogowy = 0;
graniczna_cyfrowy = 0;
dokladnosc = 1e-5;
for i=1:n
    if (graniczna_analogowy == 0 && abs( mag2db(abs(H_s(i))) ) < dokladnosc)
        graniczna_analogowy = i;
    end
    
    if (graniczna_cyfrowy == 0 && abs( mag2db(abs(H_d(i))) ) < dokladnosc)
        graniczna_cyfrowy = i;
    end
    
    if (graniczna_analogowy ~= 0 && graniczna_cyfrowy ~=0)
        break;
    end
end

disp(['Czestotliwosc graniczna dla filtru pseudo analogowego to: ', num2str(graniczna_analogowy), '[Hz]']);
disp(['Czestotliwosc graniczna dla filtru cyfrowego to: ', num2str(graniczna_cyfrowy), '[Hz]']);
disp(['Roznica pomiedzy czestotliwosciami granicznymi: ', num2str(abs(graniczna_analogowy - graniczna_cyfrowy)), '[Hz]']);


%% sygna³ cyfrowy o czasie trwania 1 s, czêstotliwoœci próbkowania fs=16 kHz, 
%  z³o¿ony z sumy dwóch harmonicznych o czêstotliwoœciach odpowiednio: 1209 i 1272 Hz

f_s = 16000;
t = 0:1/f_s:(n-1)/f_s;  %dziedzina sygnalu x [s]

f_1 = 1209; A_1 = 1; p_1 = 0; x_1 = A_1*cos(2*pi * f_1 * t + p_1);
f_2 = 1272; A_2 = 1; p_2 = 0; x_2 = A_2*cos(2*pi * f_2 * t + p_2);
x = x_1 + x_2;

figure(2);
plot(t(1, 1:n/2), x(1, 1:n/2), 'c-.');
title('Sygnal cyfrowy skladajacy sie z dwoch harmonicznych');
legend('x');
xlabel('czas [s]');


%% Cyfrowa filtracja sygnalu - implementacja wlasna

X = fft(x);                                                 % przejscie na dziedzine czestotliwosciowa
X_filtered = X .* H_d;                                      % flitracja to mnozenie w dziedzinie czestotliwosciowej
x_filtered = ifft(X(1, 1:n/2).*H_d(1, 1:n/2)) * max(x);     % powrot do dziedziny czasowej + przeskalowanie znormalizowanej osi rzednych

%% Rysuje wykres sygnalu po filtracji
figure(3)
plot(t(1, 1:n/2), abs(x_filtered), 'g-');
title('Sygnal po odfiltrowaniu jednej z dwoch harmonicznych');
legend('x po filtracji');
xlabel('czas [s]');
%UWAGA: Spadek mocy o polowe (3dB).

%% Cyfrowa filtracja sygnalu - funkcja filter
x_filtered2 = filter(b_d, a_d, x);

%% Rysuje wykres sygnalu po filtracji
figure(4)
plot(t(1, 1:n/2), real(x_filtered2(1, 1:n/2)), 'k--');
title('Sygnal po odfiltrowaniu jednej z dwoch harmonicznych');
legend('x po filtracji');
xlabel('czas [s]');