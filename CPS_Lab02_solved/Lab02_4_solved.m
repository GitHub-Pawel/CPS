close all;
clear all;
clc;

[x, fs] = audioread('mowa.wav');
%sound(x, fs);

figure(1);
hold all;
plot(x, 'k-.');
title('mowa.wav');
legend('x');
xlabel('numer probki');

N=256;
%wybieram 10 fragmentow
x_(1, :)=x(192:192+N-1);
x_(2, :)=x(1187:1187+N-1);
x_(3, :)=x(2395:2395+N-1);
x_(4, :)=x(4252:4252+N-1);
x_(5, :)=x(5000:5000+N-1);
x_(6, :)=x(6583:6583+N-1);
x_(7, :)=x(9204:9204+N-1);
x_(8, :)=x(10410:10410+N-1);
x_(9, :)=x(10980:10980+N-1);
x_(10, :)=x(12200:12200+N-1);

%obliczam macierz analizy DCT-II
s(1:N)=sqrt(2/N);   %s_0
s(1)=sqrt(1/N);     %s_k, k~=0

A=zeros(N);     
for k=1:N           %kolejne wiersze macierzy
    for n=1:N       %kolejne kolumny macierzy
        A(k,n) = s(k) * cos(pi * (k-1)/N * (n-1 + 0.5));
                    %k oraz n nale¿y przesun¹æ o 1, aby k=n=0...N-1
    end
end

%przeprowadzam analize dla kazdego z 10 fragmentow sygnalu
for i=1:10
    y(i, :)=A*x_(i, :)';
end

%skalowanie [Hz]
t = 4;            %czas trwania [s]
n = length(x);     %liczba probek [1]
Fs = n/t;    %czestotliwosc [Hz]
dziedzina = (0:(n-1))*Fs/n;

for i=1:10
    figure(2);
    ax1 = subplot(2,1,1); % top subplot
    plot(dziedzina(i:i+N-1), x_(i, :), 'k-');
    title([num2str(i), ' - sygnal']);
    xlabel('[Hz]');
    ax2 = subplot(2,1,2); % bottom subplot
    plot(dziedzina(i:i+N-1), y(i, :), 'r-');
    title([num2str(i), ' - analiza']);
    xlabel('[Hz]');
    pause;
end