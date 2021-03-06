close all;
clear all;
clc;



N = 100;
f_s = 1000;   %czestotliwosc probkowania [Hz]

f_1 = 50; f_2 = 100; f_3 = 150; %czestotliwosci sygnalow skladowych
A_1 = 50; A_2 = 100; A_3 = 150; %amplitudy sygnalow skladowych

sin_1 = A_1*sin(2*pi * f_1 * [0:1/f_s:(N-1)/f_s]);
sin_2 = A_2*sin(2*pi * f_2 * [0:1/f_s:(N-1)/f_s]);
sin_3 = A_3*sin(2*pi * f_3 * [0:1/f_s:(N-1)/f_s]);

x = sin_1 + sin_2 + sin_3;



%buduje macierze A=DCT i S=IDCT dla N=100:
s(1:N)=sqrt(2/N);   %s_0
s(1)=sqrt(1/N);     %s_k, k~=0

A=zeros(N);
for k=1:N           %kolejne wiersze macierzy
    for n=1:N       %kolejne kolumny macierzy
        A(k,n) = s(k) * cos(pi * (k-1)/N * (n-1 + 0.5));
        %k oraz n nale�y przesun�� o 1, aby k=n=0...N-1
    end
end

S = A';         %macierz odwrotna (syntezy) do A (DCT-II)
%disp(S*A);     %Sprawdzenie czy macierz jest jednostkowa



%Wy�wietlanie w p�tli warto�ci wszystkich wierszy macierzy A i kolumn macierzy S
%for i=1:N
%   for j=1:N
%       disp([A(i,j), S(j,i)]);
%   end
%end

%Wy�wietlanie przebiegow na jednym wykresie
%for i=1:N
%    figure(1);
%    hold off;
%    plot(A(i,:), 'b-o');
%    hold on;
%    plot(S(:, i)', 'r--x');
%    title('Porownanie wierszy macierzy DCT i kolumn macierzy IDCT');
%    legend('Wiersze macierzy A', 'Kolumny macierzy S');
%    text(100,0,['  i =', num2str(i)]);
%    pause;
%end



figure(1);
hold all;

y=A*x';
f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]
plot(f, y, 'c-o');

f_2 = 105;
sin_2 = A_2*sin(2*pi * f_2 * [0:1/f_s:(N-1)/f_s]);
x = sin_1 + sin_2 + sin_3;
y=A*x';
plot(f, y, 'r-*');

title('Porownanie wartosci y dla roznych czestotliwosci skladowych sygnalu x')
legend('y - dla f_2 = 100Hz', 'y - dla f_2 = 105Hz');
xlabel('[Hz]');






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

%x=rand(1,N);   %generowanie sygnalu losowego
X=A*x';         %analiza
x_s=S*X;        %rekonstrukcja (synteza)

%Sprawdzenie czy sygnaly x i x_s sa identyczne
flag = 0;
for i=1:N
    if (abs(x(i) - x_s(i)) > 1e-11)
        flag = 1;
        break;
    end
end

if flag == 1
    disp('Sygnaly x i x_s nie sa identyczne - rekonstrukcja nie jest perfekcyjna');
else
    disp('Sygnaly x i x_s sa identyczne - rekonstrukcja jest perfekcyjna');
end

figure(2);
hold all;
plot(x', 'b-o');
plot(x_s, 'r--x');
title('Porownanie sygnalu oryginalnego z sygnalem zrekonstruowanym');
legend('x - oryginalny', 'x_s - zrekonstruowany');
grid;