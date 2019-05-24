%DIT - (ang Decimation in time) - Podzia³ w dziedzinie czasu
%mag2db(y) = 20log(y);

close all;
clear all;
clc;

%==========Wygeneruj sygna³ x, losowy, o d³ugoœci 1024 próbek==========
N = 1024;           %Dlugosc sygnalu
x = rand(1, N);     %Sygnal losowy
f_s = 1000;         %czestotliwosc probkowania [Hz]

W_N = exp(1i*(2*pi/N)); 
c(1:N) = (W_N).^(-((1:N)-1));
k=0:N/2-1;


%==========Oblicz X za pomoc¹ funkcji DFT(...)==========
X = fft(x);
f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]


%==========Nastêpnie wyznacz X_fft za pomoc¹ (6): X_fft=X_1+c*X_2==========
x_1 = x(1:2:length(x));
x_2 = x(2:2:length(x));

X_1 = fft(x_1);
X_2 = fft(x_2);

%k=0:N/2-1;

X_fft(k+1) = X_1 + c(k+1).*X_2;
X_fft(k+1+N/2) = X_1 + c(k+1+N/2).*X_2;

f_fft = (0:(N/2)-1)*f_s/(N);  %przeskalowanie dziedziny [Hz]


%=========Nastêpnie widma X1 oraz X2 wyznacz ponownie za pomoc¹ (2)========
x_11 = x_1(1:2:length(x_1));
x_12 = x_1(2:2:length(x_1));
x_21 = x_2(1:2:length(x_2));
x_22 = x_2(2:2:length(x_2));

X_11 = fft(x_11); 
X_12 = fft(x_12); 
X_21 = fft(x_21); 
X_22 = fft(x_22); 

l=0:N/4-1;
%c(k+1) = (W_N).^(-k);  (-(k + N/2));

X_q1(k+1) = X_11 + c(k+1).*X_12;
X_q1(k+1 + N/2) = X_11 + c(k+1 + N/2).*X_12;

X_q2(k+1) = X_21 + c(k+1).*X_22;
X_q2(k+1 + N/2) = X_21 + c(k+1 + N/2).*X_22;

X_fft2(k+1) = X_q1 + c(k+1).*X_q2;
X_fft2(k+1+N/2) = X_q1 + c(k+1+N/2).*X_q2;

%----------Rysuje wykres----------
figure(1);
plot(f, mag2db(abs(X(1:N))), 'r-');
%figure(2);
hold all;
plot(f, mag2db(abs(X_fft2)), 'c--');

