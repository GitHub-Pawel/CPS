close all;
clear all;
clc;

%DFT - Discrete Fourier Transform
%DFT jest transformacj¹ ortogonaln¹ w przestrzeni wektorowej N-wymiarowej

%----------dane----------
N = 100;
M = 100;

f_s = 1000; %czestotliwosc probkowania [Hz]
f_1 = 125; f_2=200; %czestotliwosci skladowe [Hz]
A_1 = 100; A_2=200; %amplitudy skladowe [Hz]
p_1 = pi/7; p_2 = pi/11; %fazy skladowe [rad]
t = 0:1/f_s:(N-1)/f_s;  %dziedzina sygnalu x [s]
f=0:0.25:1000;

%----------Wyznaczam macierz A transformacji DFT----------
W_N = exp(1i*(2*pi/N));

A = zeros(100);
for k = 0:(N-1)
    for n = 0:(N-1)
        A(k+1, n+1) = (1/sqrt(N))*W_N^(-k*n);
    end
end

%----------Wyznaczam sygnal x----------
%bezposrednio:
%x = A_1*cos(2*pi * f_1 * t + p_1) + A_2*cos(2*pi * f_2 * t + p_2);

%ze wzoru: cos(a+b)=cos(a)cos(b)-sin(a)sin(b):
x_cos_1 = A_1*cos(2*pi*f_1*t)*cos(p_1);
x_cos_2 = A_2*cos(2*pi*f_2*t)*cos(p_2); 
x_sin_1 = -A_1*sin(2*pi*f_1*t)*sin(p_1);
x_sin_2 = -A_2*sin(2*pi*f_2*t)*sin(p_2);
x = x_cos_1 + x_cos_2 + x_sin_1 + x_sin_2;

%----------Wyznaczam sygnal X_1, ktore jest widmem sygnalu x----------
X_1 = A*x';            %DFT o d³ugoœci N

fx_1 = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz



%----------Wykonuje skalowanie X_2----------
x_z = zeros(1, N+M);
x_z(1, 1:N) = x;                %Zwiekszam rozdzielczosc sygnalu x

X_2 = fft(x_z)./(N+M);          %DFT z dodaniem zer

fx_2 = (0:(N+M)-1)*f_s/(N+M);   %przeskalowanie dziedziny [Hz]


%----------Wyznaczam X_3 stosujac wzor na DtFT(x)----------
X_3 = zeros(1, length(f));
for i=1:length(f)
    for n=1:N
        X_3(i) = X_3(i) + 1/N*x(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
    end
end

fx_3 = (0:length(f)-1)*f_s/length(f);  %przeskalowanie dziedziny [Hz]


%----------Rysuje wykresy widm----------
figure(1);
subplot(3, 1, 1);
hold all;
plot(fx_1, abs(X_1), 'c-'); 
plot([f_s/2 f_s/2], [min(abs(X_1)) max(abs(X_1))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_1');
legend('DFT o d³ugoœci N', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 2);
hold all;
plot(fx_2, abs(X_2), 'b-');
plot([f_s/2 f_s/2], [min(abs(X_2)) max(abs(X_2))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_2');
legend('DFT z dodaniem zer', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 3);
hold all;
plot(fx_3, abs(X_3), 'r-');
plot([f_s/2 f_s/2], [min(abs(X_3)) max(abs(X_3))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_3');
legend('DtFT', 'asymptota pionowa');
xlabel('f=0:0.25:1000 Hz');









%----------Zmiana czestotliwosci sygnalu---------
df = 0.25;
f = -2*f_s:df:2*f_s;
%f = 0:df:f_s/2;    %Wersja rownowazna powyzszej


%----------Ponownie wyznaczam X_3 stosujac wzor na DtFT(x)----------
X_3 = zeros(1, length(f));
for i=1:length(f)
    for n=1:N
        X_3(i) = X_3(i) + 1/N*x(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
    end
end

fx_3 = (0:length(f)-1)*f_s/length(f);  %przeskalowanie dziedziny [Hz]


%----------Rysuje wykresy widm----------
figure(2);
subplot(3, 1, 1);
hold all;
plot(fx_1, abs(X_1), 'c-'); 
plot([f_s/2 f_s/2], [min(abs(X_1)) max(abs(X_1))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_1');
legend('DFT o d³ugoœci N', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 2);
hold all;
plot(fx_2, abs(X_2), 'b-');
plot([f_s/2 f_s/2], [min(abs(X_2)) max(abs(X_2))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_2');
legend('DFT z dodaniem zer', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 3);
hold all;
plot(fx_3, abs(X_3), 'r-');
plot([f_s/2 f_s/2], [min(abs(X_3)) max(abs(X_3))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_3');
legend('DtFT', 'asymptota pionowa');
xlabel('f = -2000:0.25:2000 Hz');









%----------Zmiana czestotliwosci sygnalu---------
df = 0.25;
f = 0:df:f_s/2;


%----------Ponownie wyznaczam X_3 stosujac wzor na DtFT(x)----------
X_3 = zeros(1, length(f));
for i=1:length(f)
    for n=1:N
        X_3(i) = X_3(i) + 1/N*x(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
    end
end

fx_3 = (0:length(f)-1)*f_s/length(f);  %przeskalowanie dziedziny [Hz]


%----------Rysuje wykresy widm----------
figure(3);
subplot(3, 1, 1);
hold all;
plot(fx_1, abs(X_1), 'c-'); 
plot([f_s/2 f_s/2], [min(abs(X_1)) max(abs(X_1))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_1');
legend('DFT o d³ugoœci N', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 2);
hold all;
plot(fx_2, abs(X_2), 'b-');
plot([f_s/2 f_s/2], [min(abs(X_2)) max(abs(X_2))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_2');
legend('DFT z dodaniem zer', 'asymptota pionowa');
xlabel('[Hz]');

subplot(3, 1, 3);
hold all;
plot(fx_3, abs(X_3), 'r-');
%plot([f_s/2 f_s/2], [min(abs(X_3)) max(abs(X_3))], 'k:');
title('Charakterystyka amplitudowo-czestotliwosciowa widma X_3');
legend('DtFT'); %, 'asymptota pionowa');
xlabel('f = 0:0.25:500 Hz');