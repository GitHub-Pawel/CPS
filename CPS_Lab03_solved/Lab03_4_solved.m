clear all;
close all;
clc;

%----------Wczytuje plik ze sygnalami ADSL----------
load('lab_03.mat');


%----------Dane----------
M = 32;                 %Dlugosc prefiksu
N = 512;                %Dlugosc wlasciwej ramki
K = 8;                  %Liczba ramek
f_s = 2208000;          %Czêstotliwoœæ próbkowania w ADSL
f = (0:(N-1))*f_s/N;    %Skalowanie czestotliwosci [Hz]

id = mod(296892, 16)+1;
disp(['Numer wylosowanego sygnalu: ', num2str(id)]);    %13


%---------Usuwam prefiksy z sygnalu ADSL--------
x = zeros(1, K*N);
j = 0;
for k=0:K-1
    for i=1:(M+N)
        if (i>32)
            j = j+1;
            x(j) = x_13(i + (M+N)*k);
        end
    end
end


%---Wyznaczam N-punktowe DTF (FFT) kazdej ramki (po usunieciuprefiksu)---
X = zeros(K, N);    %Kolejne wiersze macierzy X, sa widmami kolejnych ramek
for k=1:K
    X(k, :) = fft(x((k-1)*N+1:N*k));
end


%Wyznaczam na podstawie wykresu, ktore harmoniczne byly uzywane w sygnale
for i=1:K
    figure(i);
    subplot(3,1,2);
    plot(f, real(X(i,:)), 'r-');
    title(['Czesc rzeczywista ramki nr: ', num2str(i), ' (ile jest cosinusa dla danej czestotliwosci)']);
    legend(['real(X(', num2str(i), ', : )']);
    xlabel('[Hz]');
    
    subplot(3,1,3);
    plot(f, imag(X(i,:)), 'b-');
    title(['Czesc urojona ramki nr ', num2str(i), ' (ile jest sinusa dla danej czestotliwosci)']);
    legend(['imag(X(', num2str(i), ', : )']);
    xlabel('[Hz]');
    
    subplot(3,1,1);
    plot(x((i-1)*N+1:N*i), 'c-');
    title(['Analizowana ramka nr ', num2str(i),]);
    legend(['x']);
end