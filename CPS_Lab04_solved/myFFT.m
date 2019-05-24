close all;
clear all;
clc;

%============================Dane==============================
N = 1024;           %Dlugosc sygnalu
f_s = 1000;         %czestotliwosc probkowania [Hz]
f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]


%===================Pobranie wektorow===================
load('x.mat', '-ascii');
load('xcpp.dat');


%===================Obliczenie transformat===================
%X = fft(x(1, :) + 1i*x(2, :));
%X_cpp = fft(xcpp(1, :) + 1i*xcpp(2, :));
X = fft(x(:, 1) + 1i*x(:, 2));
X_cpp = fft(xcpp(:, 1) + 1i*xcpp(:, 2));


%===================Porownanie widm X oraz X_cpp===================
figure(1);
hold all;
plot(f, mag2db(abs(X)), 'r-');
plot(f, mag2db(abs(X_cpp)), '.');
title('Porownanie widm sygnalu x.mat oraz xcpp.dat');
legend('X = |fft(x.mat)|', 'X_cpp = |fft(xcpp.dat)|');
xlabel('[Hz]');
ylabel('[dB]');



%----------Sprawdzenie czy sygnaly x i xcpp sa identyczne----------
flag = false;
for i=1:N
    if (abs(x(i) - xcpp(i)) > 1e-5)
        flag = true;
    end
end

if flag == true
    disp('Sygnaly x i xcpp NIE sa identyczne');
else
    disp('Sygnaly x i xcpp sa identyczne');
end




%----------Sprawdzenie czy widma X_fft1 i X_1 sa identyczne----------
flag = false;
for i=1:N
    if (abs(abs(X(i)) - abs(X_cpp(i))) > 1e-5)
        flag = true;
    end
end

if flag == true
    disp('Widma X i X_cpp NIE sa identyczne');
else
    disp('Widma X i X_cpp sa identyczne');
end