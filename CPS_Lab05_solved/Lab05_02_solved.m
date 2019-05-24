%% 2. Filtr Butterworth LP
% LP - Low Pass - dolnoprzepustowy
% Zalety filtru Butterwortha:
%   1) nie ma zafalowañ w paœmie przepustowym i zaporowym
%   2) ma charakterystykê fazow¹ najbardziej zbli¿on¹ do liniowej
%      co gwarantuje, ¿e uk³ad nie zmienia na wyjœciu kszta³tu
%      sygna³u zawartego w paœmie przepustowym.
% Wady filtru Butterwortha:
%   1) odznacza siê  najmniej stromymi zboczami pasm przejœciowych 

close all;
clear all;
clc;

%% Dane
w_3dB = 2*pi*100;   %rd/s
w = 0:0.1:3000;
n = length(w);
f = w.*pi/2;

%% Wykresy polozenia biegunow transmitancji
figure(1);
hold all;
plot(real(p(2)), imag(p(2)), 'ro');
plot(real(p(4)), imag(p(4)), 'gx');
plot(real(p(6)), imag(p(6)), 'c*');
plot(real(p(8)), imag(p(8)), 'bs');
axis([-800 800 -800 800]);
title('Polozenie biegunow transmitancji dla filtrow o kolejnych zedach');
legend('p(2)', 'p(4)', 'p(6)', 'p(8)');
xlabel('real part');
ylabel('imaginary part');


%% Wyznaczam filtry Butterwortha
H_2 = (polyval(poly(p(2)), zeros(1, n))) ./ (polyval(poly(p(2)), 1i*w));
H_4 = (polyval(poly(p(4)), zeros(1, n))) ./ (polyval(poly(p(4)), 1i*w));
H_6 = (polyval(poly(p(6)), zeros(1, n))) ./ (polyval(poly(p(6)), 1i*w));
H_8 = (polyval(poly(p(8)), zeros(1, n))) ./ (polyval(poly(p(8)), 1i*w));


%% Wykres charakterystyki amplitudowo-czestotliwosciowej liniowo
figure(2);
hold all;
plot(f, abs(H_2), 'r');    semilogx(f, abs(H_2), 'c--');
plot(f, abs(H_4), 'g');    semilogx(f, abs(H_4), 'k--');
plot(f, abs(H_6), 'c');    semilogx(f, abs(H_6), 'r--');
plot(f, abs(H_8), 'b');    semilogx(f, abs(H_8), 'g--');
title('Cha-ka a-cz kolejnych transmitancji w funkcji f wyskalowanej liniowo/logarytmicznie');
legend('|H_2(jw)|', '|H_4(jw)|', '|H_6(jw)|', '|H_8(jw)|');
xlabel('[Hz]');


%% Wykres charakterystyki amplitudowo-czestotliwosciowej logarytmicznie
figure(3);
hold all;
plot(f, mag2db(abs(H_2)), 'r');    semilogx(f, mag2db(abs(H_2)), 'c--');
plot(f, mag2db(abs(H_4)), 'g');    semilogx(f, mag2db(abs(H_4)), 'k--');
plot(f, mag2db(abs(H_6)), 'c');    semilogx(f, mag2db(abs(H_6)), 'r--');
plot(f, mag2db(abs(H_8)), 'b');    semilogx(f, mag2db(abs(H_8)), 'g--');
title('Cha-ka a-cz kolejnych transmitancji w funkcji f wyskalowanej liniowo/logarytmicznie');
legend('20log|H_2(jw)|', '20log|H_4(jw)|', '20log|H_6(jw)|', '20log|H_8(jw)|');
xlabel('[Hz]');
ylabel('[dB]');


%% Wykres charakterystyki fazowo-czestotliwosciowej liniowej
figure(4);
hold all;
plot(f, unwrap(angle(H_2)), 'r');
plot(f, unwrap(angle(H_4)), 'g');
plot(f, unwrap(angle(H_6)), 'c');
plot(f, unwrap(angle(H_8)), 'b');
title('Cha-ka f-cz kolejnych transmitancji w funkcji f wyskalowanej liniowo');
legend('unwrap(angle(H_2(jw))', 'unwrap(angle(H_4(jw))', 'unwrap(angle(H_6(jw))', 'unwrap(angle(H_8(jw))');
xlabel('[Hz]');


%% Wyznaczam odp impulsowa filtru dla N=4 oraz jego odpowiedz na skok jednostkowy
[a,b] = zp2tf([],p(2)',1); % zera, bieguny -> wspolczynniki wielomianow
H = tf(a,b); %transmitancja

%% Wykres odp impulsowej filtru dla N = 4
figure(5);
impulse(H, 'r');
legend('impulse(H)');

%% Wykres odpowiedzi filtru dla N=4 na skok jednostkowy
figure(6);
step(H, 'b');
legend('step(H)');

%% ========================================================================
