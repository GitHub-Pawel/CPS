%Zadanie 1A
clear all;
clc;

T = 0.1;            %czas trwania sygna�u [s]
A = 230;            %amplituda [V]
f = 50;             %czestotliwosc [Hz]
Usk = A*sqrt(2);    %warto�� skuteczna [V]

f_s1 = 10000;  %cz�stotliwo�� probkowanie 10kHz (pseudoanalog)
f_s2 = 500;    %cz�stotliwo�� probkowanie 500Hz
f_s3 = 200;    %cz�stotliwo�� probkowanie 200Hz

t_s1 = 0 : 1/f_s1 : T;    %wektor czas�w pr�bkowania dla sygna�u 10kHz
t_s2 = 0 : 1/f_s2 : T;    %wektor czas�w pr�bkowania dla sygna�u 500Hz
t_s3 = 0 : 1/f_s3 : T;    %wektor czas�w pr�bkowania dla sygna�u 200Hz

sin_s1 = Usk * sin(2*pi * f * t_s1);  %10kHz
sin_s2 = Usk * sin(2*pi * f * t_s2);  %500Hz
sin_s3 = Usk * sin(2*pi * f * t_s3);  %200Hz


figure(1);
hold all;

plot(t_s1, sin_s1, '-b');
plot(t_s2, sin_s2, 'o-r');
plot(t_s3, sin_s3, 'x-k');

title('Pr�bkowanie sinusa 50Hz 230V');
legend('10kHz','500Hz','200Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');
