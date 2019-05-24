%Lab04 - Zadanie 2.1 - Transformata Fouriera sygna³ów rzeczywistych
%mag2db(y) = 20log(y);

close all;
clear all;
clc;

%=======Przyjmij N=1024, wygeneruj losowe dane o rozk³adzie normalnym======
N = 1024;           %Dlugosc sygnalu
f_s = 1000;         %czestotliwosc probkowania [Hz]

k = 1:N-1;
l = 1:N/2-1;

x_1 = zeros(1, N);
x_2 = zeros(1, N);
for i=1:N
    x_1(i) = normrnd(0, 1);   %Pierwszy sygnal losowy o rozkladzie normalnym
    x_2(i) = normrnd(0, 1);   %Drugi sygnal losowy o rozkladzie normalnym
end


%==========Wykonaj transformatê zespolon¹ dla porównania==========
tic;

X_fft1 = fft(x_1);
X_fft2 = fft(x_2);

t1 = toc;
disp(['Czas t_1 obliczen widm za pomoca FFT(...): ', num2str(t1), 's']);

f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]


%==================================================================
tic;

y = x_1 + 1i*x_2;

Y = fft(y);
Y_r = real(Y);
Y_i = imag(Y);

X_1r(1) = Y_r(1);
X_2r(1) = Y_i(1);
X_1i(1) = 0;
X_2i(1) = 0;

X_1r(k+1) = 0.5*(Y_r(k+1) + Y_r(N-k+1));
X_2r(k+1) = 0.5*(Y_i(k+1) + Y_i(N-k+1));

X_1i(k+1) = 0.5*(Y_i(k+1) - Y_i(N-k+1));
X_2i(k+1) = 0.5*(Y_r(N-k+1) - Y_r(k+1));

X_1 = X_1r + 1i*X_1i;
X_2 = X_2r + 1i*X_2i;

t2 = toc;
disp(['Czas obliczen t_2 widm z wykorzystaniem symetrii: ', num2str(t2), 's']);
disp('----------------------------------------------------------');
disp(['Roznica czasow t_1 - t_2 = ', num2str(t1-t2), 's']);

%------------------Porownuje sygnaly na wykresach--------------------------
figure(1);
hold all;
plot(f, mag2db(abs(X_fft1)), 'r-');
plot(f, mag2db(abs(X_1)), 'c--');
title('Porownanie widm sygnalu x_1 o rozgladzie normalnym obliczonych dwoma metodami');
legend('X_1 = |fft(x_1)|', 'X_1 = |X_1_r + jX_1_i|');
xlabel('[Hz]');
ylabel('[dB]');


figure(2);
hold all;
plot(f, mag2db(abs(X_fft2)), 'b-');
plot(f, mag2db(abs(X_2)), 'g--');
title('Porownanie widm sygnalu x_2 o rozgladzie normalnym obliczonych dwoma metodami');
legend('X_2 = |fft(x_2)|', 'X_2 = |X_2_r + jX_2_i|');
xlabel('[Hz]');
ylabel('[dB]');



disp('----------------------------------------------------------');

%----------Sprawdzenie czy widma X_fft1 i X_1 sa identyczne----------
flag = false;
for i=1:N
    if (abs(abs(X_fft1(i)) - abs(X_1(i))) > 1e-10)
        flag = true;
    end
end

if flag == true
    disp('Widma X_fft1 i X_1 NIE sa identyczne');
else
    disp('Widma X_fft1 i X_1 sa identyczne');
end


%----------Sprawdzenie czy widma X_fft2 i X_2 sa identyczne----------
flag = false;
for i=1:N
    if (abs(abs(X_fft2(i)) - abs(X_2(i))) > 1e-10)
        flag = true;
    end
end

if flag == true
    disp('Widma X_fft2 i X_2 NIE sa identyczne');
else
    disp('Widma X_fft2 i X_2 sa identyczne');
end