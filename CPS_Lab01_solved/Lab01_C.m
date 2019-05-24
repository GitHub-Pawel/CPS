%Zadanie 1B sinus
clear all;
clc;

T = 1;            %czas trwania sygna³u [s]
A = 230;          %amplituda [V]
f_s = 100;        %czestotliwosc [Hz]
Usk = A*sqrt(2);  %wartoœæ skuteczna [V]

t = 0 : 1/f_s : T;  %wektor czasów próbkowania

%for f=0:5:300
%    plot(Usk*sin(2*pi * f * t));
%    title(['Próbkowanie sinusa',' ',num2str(f),'Hz, iteracja ',num2str(f/5 + 1)]);
%    legend('100 Hz');
%    xlabel('Czas [s]');
%    ylabel('Amplituda [V]');
%    pause;
%end

figure(1);
hold all;
plot(t, Usk*sin(2*pi * 5 * t), 'g-');
plot(t, Usk*sin(2*pi * 105 * t), 'r--o');
plot(t, Usk*sin(2*pi * 205 * t), 'k--x');

title('Próbkowanie sinusów 230V czêstotliwoœæ próbkowania 100Hz');
legend('Czêstotliwoœæ sinusa 5Hz','Czêstotliwoœæ sinusa 105Hz','Czêstotliwoœæ sinusa 205Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(2);
hold all;
plot(t, Usk*sin(2*pi *95 * t), 'g-');
plot(t, Usk*sin(2*pi * 195 * t), 'r--o');
plot(t, Usk*sin(2*pi * 295 * t), 'k--x');

title('Próbkowanie sinusów 230V czêstotliwoœæ próbkowania 100Hz');
legend('Czêstotliwoœæ sinusa 95Hz','Czêstotliwoœæ sinusa 195Hz','Czêstotliwoœæ sinusa 295Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(3);
hold all;
plot(t, Usk*sin(2*pi * 95 * t), 'g');
plot(t, Usk*sin(2*pi * 105 * t), 'r');

title('Próbkowanie sinusów 230V czêstotliwoœæ próbkowania 100Hz');
legend('Czêstotliwoœæ sinusa 95Hz','Czêstotliwoœæ sinusa 105Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');