%Zadanie 1B
clear all;
clc;

T = 1;            %czas trwania sygna³u [s]
A = 230;            %amplituda [V]
f = 50;             %czestotliwosc [Hz]
Usk = A*sqrt(2);    %wartoœæ skuteczna [V]

f_s1 = 10000;  %czêstotliwoœæ probkowanie 10kHz (pseudoanalog)
f_s2 = 51;     %czêstotliwoœæ probkowanie 51Hz
f_s3 = 50;     %czêstotliwoœæ probkowanie 50Hz
f_s4 = 49;     %czêstotliwoœæ probkowanie 49Hz

t_s1 = 0 : 1/f_s1 : T;    %wektor czasów próbkowania dla sygna³u 10kHz
t_s2 = 0 : 1/f_s2 : T;    %wektor czasów próbkowania dla sygna³u 51Hz
t_s3 = 0 : 1/f_s3 : T;    %wektor czasów próbkowania dla sygna³u 50Hz
t_s4 = 0 : 1/f_s4 : T;    %wektor czasów próbkowania dla sygna³u 49Hz

sin_s1 = Usk * sin(2*pi * f * t_s1);  %10kHz
sin_s2 = Usk * sin(2*pi * f * t_s2);  %51Hz
sin_s3 = Usk * sin(2*pi * f * t_s3);  %50Hz
sin_s4 = Usk * sin(2*pi * f * t_s4);  %49Hz


figure(1);
hold all;

plot(t_s1, sin_s1, 'b-');
plot(t_s2, sin_s2, 'g-o');
plot(t_s3, sin_s3, 'r-o');
plot(t_s4, sin_s4, 'k-o');

title('Próbkowanie sinusa 50Hz 230V');
legend('10kHz','51Hz','50Hz','49Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');

%==========================================================================

f_s5 = 26;     %czêstotliwoœæ probkowanie 26Hz
f_s6 = 25;     %czêstotliwoœæ probkowanie 25Hz
f_s7 = 24;     %czêstotliwoœæ probkowanie 24Hz

t_s5 = 0 : 1/f_s5 : T;    %wektor czasów próbkowania dla sygna³u 26Hz
t_s6 = 0 : 1/f_s6 : T;    %wektor czasów próbkowania dla sygna³u 25Hz
t_s7 = 0 : 1/f_s7 : T;    %wektor czasów próbkowania dla sygna³u 24Hz

sin_s5 = Usk * sin(2*pi * f * t_s5);  %26Hz
sin_s6 = Usk * sin(2*pi * f * t_s6);  %25Hz
sin_s7 = Usk * sin(2*pi * f * t_s7);  %24Hz


figure(2);
hold all;

plot(t_s1, sin_s1, 'b-');
plot(t_s5, sin_s5, 'g-o');
plot(t_s6, sin_s6, 'r-o');
plot(t_s7, sin_s7, 'k-o');

title('Próbkowanie sinusa 50Hz 230V');
legend('10kHz','26Hz','25Hz','24Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');