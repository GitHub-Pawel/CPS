%Zadanie 1A
clear all;
clc;

T = 0.1;            %czas trwania sygna³u [s]
A = 230;            %amplituda [V]
f = 50;             %czestotliwosc [Hz]
Usk = A*sqrt(2);    %wartoœæ skuteczna [V]

f_s1 = 10000;  %czêstotliwoœæ probkowanie 10kHz (pseudoanalog)
f_s2 = 500;    %czêstotliwoœæ probkowanie 500Hz
f_s3 = 200;    %czêstotliwoœæ probkowanie 200Hz

t_s1 = 0 : 1/f_s1 : T;    %wektor czasów próbkowania dla sygna³u 10kHz
t_s2 = 0 : 1/f_s2 : T;    %wektor czasów próbkowania dla sygna³u 500Hz
t_s3 = 0 : 1/f_s3 : T;    %wektor czasów próbkowania dla sygna³u 200Hz

sin_s1 = Usk * sin(2*pi * f * t_s1);  %Sygna³ spróbkowany 10kHz
sin_s2 = Usk * sin(2*pi * f * t_s2);  %Sygna³ spróbkowany 500Hz
sin_s3 = Usk * sin(2*pi * f * t_s3);  %Sygna³ spróbkowany 200Hz

Sa_s1 = sinc(2*pi * f * t_s1);  %Sampling do odtworzenia sygna³u spróbkowanego 10kHz
Sa_s2 = sinc(2*pi * f * t_s2);  %Sampling do odtworzenia sygna³u spróbkowanego 500Hz
Sa_s3 = sinc(2*pi * f * t_s3);  %Sampling do odtworzenia sygna³u spróbkowanego 200Hz

%Czêœæ kodu z instrukcji
xhat = zeros(1, length(t_s2));
for idx = 1:length(t_s2)
        for n = 1:length(t_s2)
            xhat(idx) = xhat(idx) + sin_s2(n)*sinc(pi*(idx-n));     
        end
end

figure(1);
hold all;

plot(t_s2, sin_s2, 'r-x');
plot(t_s2, xhat, 'g--o');
plot(t_s2, Usk*Sa_s2);

title('Odtwarzanie sinusa 50Hz 230V');
legend('Sygna³ spróbkowany 500Hz', 'Sygna³ odtworzony', 'Sampling wykorzystany do odtworzenia');
xlabel('Czas [s]');
ylabel('Amplituda [V]');

%==========================================================================

%Czêœæ kodu z instrukcji
xhat2 = zeros(1, length(t_s3));
for idx = 1:length(t_s3)
        for n = 1:length(t_s3)
            xhat2(idx) = xhat2(idx) + sin_s3(n)*sinc(pi*(idx-n));     
        end
end

figure(2);
hold all;

plot(t_s3, sin_s3, 'r-x');
plot(t_s3, xhat2, 'g--o');
plot(t_s3, Usk*Sa_s3);

title('Odtwarzanie sinusa 50Hz 230V');
legend('Sygna³ spróbkowany 200Hz', 'Sygna³ odtworzony', 'Sampling wykorzystany do odtworzenia');
xlabel('Czas [s]');
ylabel('Amplituda [V]');