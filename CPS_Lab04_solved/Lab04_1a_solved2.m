%Lab04 - Zadanie 1a - FFT sk³adane od drugiego poziomu „motylków”
%DIT - (ang Decimation in time) - Podzia³ w dziedzinie czasu
%mag2db(y) = 20log(y);

close all;
clear all;
clc;

%==========Wygeneruj sygna³ x, losowy, o d³ugoœci 1024 próbek==========
N = 1024;           %Dlugosc sygnalu
x = rand(1, N);     %Sygnal losowy
f_s = 1000;         %czestotliwosc probkowania [Hz]


%==========Oblicz X za pomoc¹ funkcji DFT(...)==========
tic;

X = fft(x);

t = toc;
disp(['Czas obliczen X za pomoca DFT(...): ', num2str(t)]);

f = (0:N-1)*f_s/N;  %przeskalowanie dziedziny [Hz]


%==========Nastêpnie wyznacz X_fft za pomoc¹ (6): X_fft=X_1+c*X_2==========
tic; 

x_1 = x(1:2:length(x));         %Probki parzyste sygnalu x
x_2 = x(2:2:length(x));         %Probki nieparzyste sygnalu x

X_1 = fft(x_1);         %Widmo DFT dla probek parzystych
X_2 = fft(x_2);         %Widmo DFT dla probek nieparzystych

k=0:N/2-1;
c(1:N) = (W(N)).^(-((1:N)-1));      %Korekta dla N probek

X_fft(k+1) = X_1 + c(k+1).*X_2;             %"Pierwsza polowa" widma X
X_fft(k+1+N/2) = X_1 + c(k+1+N/2).*X_2;     %"Druga polowa" widma X

t = toc;
disp(['Czas obliczen X = X_1 + cX_2: ', num2str(t)]);

f_fft = (0:(N/2)-1)*f_s/(N);  %przeskalowanie dziedziny [Hz]


%=========Nastêpnie widma X1 oraz X2 wyznacz ponownie za pomoc¹ (2)========
tic;

x_11 = x_1(1:2:length(x_1));    %Probki parzyste sygnalu x_1
x_12 = x_1(2:2:length(x_1));    %Probki nieparzyste sygnalu x_1
x_21 = x_2(1:2:length(x_2));    %Probki parzyste sygnalu x_2
x_22 = x_2(2:2:length(x_2));    %Probki nieparzyste sygnalu x_1

X_11 = fft(x_11);   %Widmo DFT dla probek parzystych sygnalu x_1
X_12 = fft(x_12);   %Widmo DFT dla probek nieparzystych sygnalu x_1
X_21 = fft(x_21);   %Widmo DFT dla probek parzystych sygnalu x_2
X_22 = fft(x_22);   %Widmo DFT dla probek nieparzystych sygnalu x_2

l=0:N/4-1;
c2(1:N/2) = (W(N/2)).^(-((1:N/2)-1));   %Korekta dla N/2 probek

X_q1(l+1) = X_11 + c2(l+1).*X_12;       
X_q1(l+1 + N/4) = X_11 + c2(l+1 + N/4).*X_12;

X_q2(l+1) = X_21 + c2(l+1).*X_22;
X_q2(l+1 + N/4) = X_21 + c2(l+1 + N/4).*X_22;

X_fft2(k+1) = X_q1 + c(k+1).*X_q2;          %"Pierwsza polowa" widma X
X_fft2(k+1+N/2) = X_q1 + c(k+1+N/2).*X_q2;  %"Druga polowa" widma X

t = toc;
disp(['Czas obliczen X = (X_1_1 + cX_1_2) + c(X_2_1 + cX_2_2): ', num2str(t)]);


%===Porównaj czy wynik uzyskany we wszystkich 3 sposobach jest taki sam===
figure(1);
hold all;
plot(f, mag2db(abs(X)), 'r-');
plot(f, mag2db(abs(X_fft)), 'c--');
plot(f, mag2db(abs(X_fft2)), 'k.');
title("Porownanie wynikow uzyskanych we wszystkich 3 sposobach");
legend('X = DFT(...)', 'X_f_f_t = X_1 + cX_2', 'X_f_f_t = (X_1_1 + cX_1_2) + c(X_2_1 + cX_2_2)');
xlabel('[Hz]');
ylabel('[dB]');




%----------Sprawdzenie czy widma X i X_fft sa identyczne----------
flag = false;
for i=1:N
    if (abs(abs(X(i)) - abs(X_fft(i))) > 1e-10)
        flag = true;
    end
end

if flag == true
    disp('Widma X i X_fft NIE sa identyczne');
else
    disp('Widma X i X_fft sa identyczne');
end



%----------Sprawdzenie czy widma X i X_fft2 sa identyczne----------
flag = false;
for i=1:N
    if (abs(abs(X(i)) - abs(X_fft2(i))) > 1e-10)
        flag = true;
    end
end

if flag == true
    disp('Widma X i X_fft2 NIE sa identyczne');
else
    disp('Widma X i X_fft2 sa identyczne');
end