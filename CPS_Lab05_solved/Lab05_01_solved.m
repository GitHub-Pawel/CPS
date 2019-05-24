%% 1. Projektowanie metod¹ zer i biegunów 
%LTI - Linear Time-Invariant - Liniowy niezmienny w czasie

close all;
clear all
clc;

%% ============Dane============
% p - wektor pierwiastkow wielomianu z mianownika
p(1) = -0.5 + 1i*9.5;
p(3) = -1 + 1i*10;
p(5) = -0.5 + 1i*10.5 ;
for i=2:2:length(p)+1
    p(i) = p(i-1)';
end
N = length(p);

% z - wektor pierwiastkow wielomianu z licznika
z(1) = 5i;
z(3) = 15i;
for i=2:2:length(z)+1
    z(i) = z(i-1)';
end
M = length(z);

a_N = 1;
b_M = 1;

n = 20;
f_s = 1000;
w = 0:1/f_s:n-1/f_s;

%% ============Obliczanie transmitancji ze wzoru (1): SPOSOB 1============
for i=1:N
    a_N = a_N .* (1i*w - p(i));
end

for i=1:M
    b_M = b_M .* (1i*w - z(i));
end

H = b_M ./ a_N;

%% ============Obliczanie transmitancji ze wzoru (1): SPOSOB 2============
% p = poly(roots)  - dla zadanego wektora pierwiastkow wielomianu 'roots', 
%                    zwraca wektor jego wspolczynnikow.
% y = polyval(p,x) - dla zadanego wektora wspolczynnikow oraz wektora
%                    argumentow, zwraca sprobkowany wielomian.

a_N = 1;
b_M = 1;

H = (a_N*polyval(poly(z), 1i*w)) ./ (b_M*polyval(poly(p), 1i*w));

%% ============wykres zer i biegunow na p³aszczyznie zespolonej============
figure(1);
hold all;
plot(real(z),imag(z), 'ro')
plot(real(p),imag(p),'b*');
legend('Zera','Bieguny');
xlabel('Real');
ylabel('Imag');
grid;
title('Zera i bieguny');


%% ============Wykres cha-ki a-cz w skali liniowej============
figure(2);
plot(w,abs(H), 'r');
legend('y = |H(jw)|');
xlabel('x');
ylabel('y');
title('Charakterystyka amplitudowo-czestotliwosciowa H(jw) w skali liniowej');


%% ============Wykres cha-ki a-cz w skali logarytmicznej============
figure(3);
plot(w,mag2db(abs(H)), 'b');
legend('y = 20*log(|H(jw)|)');
xlabel('x');
ylabel('y [Db]');
title('Charakterystyka amplitudowo-czestotliwosciowa H(jw) w skali logarytmicznej');

%% Wyznaczam maksymalne i minimalne tlumienie H(jw) w pasmie zaporowym

[max_left, id_max_left] = max(mag2db(abs(H(1 : imag( z(1) )*f_s ))));
[min_left, id_min_left] = min(mag2db(abs(H(1 : imag( z(1) )*f_s ))));

[max_right, id_max_right] = max(mag2db(abs(H(imag( z(3) )*f_s + 2: n*f_s ))));
[min_right, id_min_right] = min(mag2db(abs(H(imag( z(3) )*f_s + 2: n*f_s ))));

disp(['Max tlumienie w lewej czesci pasma zaporowego: ', num2str(max_left), '[dB] przyjmowane jest dla x = ', num2str(id_max_left/f_s)]);
disp(['Min tlumienie w lewej czesci pasma zaporowego: ', num2str(min_left), '[dB] przyjmowane jest dla x = ', num2str(id_min_left/f_s)]);

disp(['Max tlumienie w prawej czesci pasma zaporowego: ', num2str(max_right), '[dB] przyjmowane jest dla x = ', num2str(id_max_right/f_s + imag(z(3)))]);
disp(['Min tlumienie w prawej czesci pasma zaporowego: ', num2str(min_right), '[dB] przyjmowane jest dla x = ', num2str(id_min_right/f_s + imag(z(3)))]);

%% Normalizacja transmitancji wzmocnienia w pasmie podstawowym do 1
srodek = abs(imag(z(1)) - imag(z(3)))*f_s;
disp(['Wzmmocnienie w pasmie przepustowym przed normalizacja: ', num2str(abs(H(srodek)))]);

H = H/abs(H(srodek));
disp(['Wzmmocnienie w pasmie przepustowym po normalizacji: ', num2str(abs(H(srodek)))]);

%% ============Wykres cha-ki a-cz w skali liniowej PO NORMALIZACJI============
figure(4);
plot(w,abs(H), 'r');
legend('y = |H(jw)|');
xlabel('x');
ylabel('y');
title('Charakterystyka a-cz H(jw) w skali liniowej po normalizacji');
grid;

%% ========Wykres cha-ki a-cz w skali logarytmicznej PO NORMALIZACJI=======
figure(5);
plot(w,mag2db(abs(H)), 'b');
legend('y = 20*log(|H(jw)|)');
xlabel('x');
ylabel('y [Db]');
title('Charakterystyka a-cz H(jw) w skali logarytmicznej po normalizacji');
grid;

%% ============Wykres cha-ki f-cz w skali liniowej PO NORMALIZACJI=========
figure(6);
plot(w, unwrap(angle(H)), 'k');
legend('y = unwrap(angle(H(jw)))');
xlabel('x');
ylabel('y');
title('Charakterystyka fazowo-czestotliwosciowa H(jw) w skali liniowej po normalizacji');
grid;

%% ========================================================================