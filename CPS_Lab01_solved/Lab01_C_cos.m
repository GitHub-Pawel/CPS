%Zadanie 1B cosinus
clear all;
clc;

T = 1;            %czas trwania sygnału [s]
A = 230;          %amplituda [V]
f_s = 100;        %czestotliwosc [Hz]
Usk = A*sqrt(2);  %wartość skuteczna [V]

t = 0 : 1/f_s : T;  %wektor czasów próbkowania

%for f=0:5:300
%    plot(Usk*cos(2*pi * f * t));
%    title(['Próbkowanie cosinusa',' ',num2str(f),'Hz, iteracja ',num2str(f/5 + 1)]);
%    legend('100 Hz');
%    xlabel('Czas [s]');
%    ylabel('Amplituda [V]');
%    pause;
%end

figure(1);
hold all;
plot(t, Usk*cos(2*pi * 5 * t), 'g-');
plot(t, Usk*cos(2*pi * 105 * t), 'r--o');
plot(t, Usk*cos(2*pi * 205 * t), 'k--x');

title('Próbkowanie cosinusów 230V częstotliwość próbkowania 100Hz');
legend('Częstotliwość cosinusa 5Hz','Częstotliwość cosinusa 105Hz','Częstotliwość cosinusa 205Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(2);
hold all;
plot(t, Usk*cos(2*pi *95 * t), 'g-');
plot(t, Usk*cos(2*pi * 195 * t), 'r--o');
plot(t, Usk*cos(2*pi * 295 * t), 'k--x');

title('Próbkowanie cosinusów 230V częstotliwość próbkowania 100Hz');
legend('Częstotliwość cosinusa 95Hz','Częstotliwość cosinusa 195Hz','Częstotliwość cosinusa 295Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(3);
hold all;
plot(t, Usk*cos(2*pi * 95 * t), 'g-*');
plot(t, Usk*cos(2*pi * 105 * t), 'r-o');

title('Próbkowanie cosinusów 230V częstotliwość próbkowania 100Hz');
legend('Częstotliwość cosinusa 95Hz','Częstotliwość cosinusa 105Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');