%Zadanie 1B sinus
clear all;
clc;

T = 1;            %czas trwania sygna�u [s]
A = 230;          %amplituda [V]
f_s = 100;        %czestotliwosc [Hz]
Usk = A*sqrt(2);  %warto�� skuteczna [V]

t = 0 : 1/f_s : T;  %wektor czas�w pr�bkowania

%for f=0:5:300
%    plot(Usk*sin(2*pi * f * t));
%    title(['Pr�bkowanie sinusa',' ',num2str(f),'Hz, iteracja ',num2str(f/5 + 1)]);
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

title('Pr�bkowanie sinus�w 230V cz�stotliwo�� pr�bkowania 100Hz');
legend('Cz�stotliwo�� sinusa 5Hz','Cz�stotliwo�� sinusa 105Hz','Cz�stotliwo�� sinusa 205Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(2);
hold all;
plot(t, Usk*sin(2*pi *95 * t), 'g-');
plot(t, Usk*sin(2*pi * 195 * t), 'r--o');
plot(t, Usk*sin(2*pi * 295 * t), 'k--x');

title('Pr�bkowanie sinus�w 230V cz�stotliwo�� pr�bkowania 100Hz');
legend('Cz�stotliwo�� sinusa 95Hz','Cz�stotliwo�� sinusa 195Hz','Cz�stotliwo�� sinusa 295Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');



figure(3);
hold all;
plot(t, Usk*sin(2*pi * 95 * t), 'g');
plot(t, Usk*sin(2*pi * 105 * t), 'r');

title('Pr�bkowanie sinus�w 230V cz�stotliwo�� pr�bkowania 100Hz');
legend('Cz�stotliwo�� sinusa 95Hz','Cz�stotliwo�� sinusa 105Hz');
xlabel('Czas [s]');
ylabel('Amplituda [V]');