close all;
clear all;
clc;

N=20;
s(1:N)=sqrt(2/N);   %s_0
s(1)=sqrt(1/N);     %s_k, k~=0

A=zeros(N);     
for k=1:N           %kolejne wiersze macierzy
    for n=1:N       %kolejne kolumny macierzy
        A(k,n) = s(k) * cos(pi * (k-1 + 0.25)/N * (n-1 + 0.5));     %podstawienie k+0.25 za k
                    %k oraz n nale¿y przesun¹æ o 1, aby k=n=0...N-1
    end
end

%==========================================================================
%Sprawdzam ortogonalnosc/ortonormalnosc macierzy A
ortogonalna = 1;
ortonormalna = 1;
for i=1:N
    for j=i:N
        if (i~=j &&  abs(sum(A(i,:).*conj(A(j,:)))) > 1e-12)
            ortogonalna = 0;
        elseif (i==j &&  abs(sum(A(i,:).*conj(A(j,:))) - 1) > 1e-12)
            ortonormalna = 0;
        end
    end
end

if (ortogonalna == 1 && ortonormalna == 1)
    disp('Macierz A jest ortonormalna');
elseif (ortogonalna == 1)
    disp('Macierz A jest ortogonalna');
else
    disp('Macierz A nie jest ortogonalna');
end

S = A';         %macierz odwrotna (syntezy) do A (DCT-II)
%disp(S*A);     %Sprawdzenie czy macierz jest jednostkowa

x=rand(1,N);   %generowanie sygnalu losowego
X=A*x';         %analiza
x_s=S*X;        %rekonstrukcja (synteza)

%Sprawdzenie czy sygnaly x i x_s sa identyczne
flag = 0;
for i=1:N
    if (abs(x(i) - x_s(i)) > 1e-12)
        flag = 1;
        break;
    end
end

if flag == 1
    disp('Sygnaly x i x_s nie sa identyczne - rekonstrukcja nie jest perfekcyjna');
else
    disp('Sygnaly x i x_s sa identyczne - rekonstrukcja jest perfekcyjna');
end

figure(1);
hold all;
plot(x', 'b-o');
plot(x_s, 'r--x');
title('Porownanie sygnalu oryginalnego z sygnalem zrekonstruowanym');
legend('x - oryginalny', 'x_s - zrekonstruowany');
grid;
