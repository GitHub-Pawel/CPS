close all;
clear all;
clc;

%DFT - Discrete Fourier Transform
%DFT jest transformacj¹ ortogonaln¹ w przestrzeni wektorowej N-wymiarowej

%----------dane----------
N = 100;

f_s = 1000; %czestotliwosc probkowania [Hz]
f_1 = 100; f_2=200; %czestotliwosci skladowe [Hz]
A_1 = 100; A_2=200; %amplitudy skladowe [Hz]
p_1 = pi/7; p_2 = pi/11; %fazy skladowe [rad]
t = 0:1/f_s:(N-1)/f_s;  %dziedzina sygnalu x [s]
f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]

f_1 = 125; %Zmiana czestotliwosci sygnalu skladowego [Hz]

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

%----------Wyznaczam sygnal X, ktore jest widmem sygnalu x----------
X = A*x';

%----------Rysuje wykres widma X----------
figure(1);

subplot(4,1,1);
plot(f, real(X), 'c-o');
hold on;
plot([f_s/2 f_s/2], [min(real(X)) max(real(X))], 'k:');
hold off;
title('Czesc rzeczywista widma X');
legend('real(X)');
xlabel('[Hz]');

subplot(4,1,2);
plot(f, imag(X), 'r-o');
hold on;
plot([f_s/2 f_s/2], [min(imag(X)) max(imag(X))], 'k:');
hold off;
title('Czesc urojona widma X');
legend('imag(X)');
xlabel('[Hz]');

subplot(4,1,3);
plot(f, abs(X), 'g-o');
hold on;
plot([f_s/2 f_s/2], [min(abs(X)) max(abs(X))], 'k:');
hold off;
title('Charakterystyka amplitudowo-czestotliwosciowa widma X');
legend('abs(X)');
xlabel('[Hz]');

subplot(4,1,4);
plot(f, unwrap(angle(X)), 'b-o');
hold on;
plot([f_s/2 f_s/2], [min(unwrap(angle(X))) max(unwrap(angle(X)))], 'k:');
hold off;
title('Charakterystyka fazowo-czestotliwosciowa widma X');
legend('unwrap(angle(X))');
xlabel('[Hz]');


%----------Sprawdzam ortogonalnosc/ortonormalnosc macierzy A----------
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


%----------Wyznaczam macierz rekonstrukcji B----------
B = A';


%----------Rekonstrukacja sygnalu x z jego widma X----------
x_r = real(B*X);


%----------Sprawdzenie czy sygnaly x i x_r sa identyczne----------
flag = 0;
for i=1:N
    if (abs(x(i) - x_r(i)) > 1e-9)
        flag = 1;
    end
end

if flag == 1
    disp('Sygnaly x i x_r nie sa identyczne - rekonstrukcja nie jest perfekcyjna');
else
    disp('Sygnaly x i x_r sa identyczne - rekonstrukcja jest perfekcyjna');
end

figure(2);
hold all;
plot(f, x', 'b-o');
plot(f, x_r, 'r--x');
title('Porownanie sygnalu oryginalnego z sygnalem zrekonstruowanym');
legend('x - oryginalny', 'x_r - zrekonstruowany');
grid;





%--------Wykorzystuje do transsformacji wbudowane funkcje MATLABA--------
X_2 = fft(x);
x_r2 = ifft(X_2);


%----------Sprawdzenie czy sygnaly x i x_r2 sa identyczne----------
flag = 0;
for i=1:N
    if (abs(x(i) - x_r2(i)) > 1e-9)
        flag = 1;
    end
end

if flag == 1
    disp('Sygnaly x i x_r2 nie sa identyczne - rekonstrukcja nie jest perfekcyjna');
else
    disp('Sygnaly x i x_r2 sa identyczne - rekonstrukcja jest perfekcyjna');
end

figure(3);
hold all;
plot(f, x', 'b-o');
plot(f, x_r2, 'r--x');
title('Porownanie sygnalu oryginalnego z sygnalem zrekonstruowanym');
legend('x - oryginalny', 'x_r2 - zrekonstruowany');
grid;




%-----------Porownuje widma wyznaczone dwoma sposobami----------
flag = 0;
for i=1:N
    if (abs(X(i) - X_2(i)/sqrt(N)) > 1e-8)
        flag = 1;
    end
end

if flag == 1
    disp('Widma X i X_2/sqrt(N) nie sa identyczne');
else
    disp('Widma X i X_2/sqrt(N) sa identyczne');
end

figure(4);
hold all;
plot(f, real(X), 'b-o');
plot(f, real(X_2/sqrt(N)), 'r--x');
title('Porownanie widm sygnalu x wyznaczonuch dwoma sposobami');
legend('X = real(A*x)', 'X_2 = real(fft(x)/sqrt(N))');
grid;
