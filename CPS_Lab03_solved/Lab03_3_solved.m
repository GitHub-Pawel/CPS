close all;
clear all;
clc;


%----------dane----------
N = 500;

f=0:0.1:500 ; %czestotliwosc sygnalu x
f_s = 1000; %czestotliwosc probkowania [Hz]
f_1 = 100; f_2=125; %czestotliwosci skladowe [Hz]
A_1 = 1; A_2=0.0001; %amplitudy skladowe [Hz]
p_1 = pi/7; p_2 = pi/11; %fazy skladowe [rad]
t = 0:1/f_s:(N-1)/f_s;  %dziedzina sygnalu x [s]



%----------Wyznaczam sygnal x----------
%x = A_1*cos(2*pi * f_1 * t + p_1) + A_2*cos(2*pi * f_2 * t + p_2);
%ze wzoru: cos(a+b)=cos(a)cos(b)-sin(a)sin(b):
x_cos_1 = A_1*cos(2*pi*f_1*t)*cos(p_1);
x_cos_2 = A_2*cos(2*pi*f_2*t)*cos(p_2); 
x_sin_1 = -A_1*sin(2*pi*f_1*t)*sin(p_1);
x_sin_2 = -A_2*sin(2*pi*f_2*t)*sin(p_2); 
%x = x_cos_1 + x_cos_2;                     %pierwsza skladowa
%x = x_sin_1 + x_sin_2;                     %druga skladowa
x = x_cos_1 + x_cos_2 + x_sin_1 + x_sin_2;  %obie skladowe


%----------Wyznaczam widmo DtFT sygnalu x-------------------
X = zeros(1, length(f));
for i=1:length(f)
    for n=1:N
        X(i) = X(i) + 1/N*x(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
    end
end

fx = (0:length(f)-1)*f_s/length(f);  %przeskalowanie dziedziny [Hz]


%----------Rysuje wykres widma X DtFT----------
figure(1);
hold all;
subplot(4,1,1);
plot(x, 'r-');
title('Analizowany wygnal x');
legend('x');

subplot(4,1,2);
plot(fx, real(X), 'b-');
title('Czesc rzeczywista sugnalu X');
legend('real(X)');
xlabel('f = 0:0.1:500 Hz');

subplot(4,1,3);
plot(fx, imag(X), 'c-');
title('Czesc urojona sugnalu X');
legend('imag(X)');
xlabel('f = 0:0.1:500 Hz');

subplot(4,1,4);
plot(fx, abs(X), 'g-');
title('Charakterystyka amplitudowo-czestotliwosciowa sugnalu X');
legend('abs(X)');
xlabel('f = 0:0.1:500 Hz');

%----------Wymnazam probki sygnalu x przez kolejne okna------------
x_r = x .* rectpuls(t, (N-1)/f_s);  %okno prostokatne
x_h = x .* hamming(N)';             %okno Hamminga 
x_b = x .* blackman(N)';            %okno Blackman
x_c100 = x .* chebwin(N, 100)';     %okno Czebyszewa (100 dB)
x_c120 = x .* chebwin(N, 120)';     %okno Czebyszewa (120 dB)



%----------Wyznaczam widma DtFT kolejnych sygnalow--------
X_r = zeros(1, length(f));
X_h = zeros(1, length(f));
X_b = zeros(1, length(f));
X_c100 = zeros(1, length(f));
X_c120 = zeros(1, length(f));
for i=1:length(f)
    for n=1:N
        X_r(i) = X_r(i) + 1/N*x_r(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
        X_h(i) = X_h(i) + 1/N*x_h(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
        X_b(i) = X_b(i) + 1/N*x_b(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
        X_c100(i) = X_c100(i) + 1/N*x_c100(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
        X_c120(i) = X_c120(i) + 1/N*x_c120(n)*exp(-1i*2*pi*f(i)/f_s*(n-1));
    end
end


%----------Rysuje wykresy modulow widm DtFT na jednym wykresie----------
figure(2);
hold all;
plot(fx, abs(X_r), 'r-');
plot(fx, abs(X_h), 'g-');
plot(fx, abs(X_b), 'c-');
plot(fx, abs(X_c100), 'b-');
plot(fx, abs(X_c120), 'k-');
title('Charakterystki amplitudowo-czestotliwosciowe przefiltrowanych sygnalow');
legend('abs(X_r_e_c_t_p_u_l_s)', 'abs(X_H_a_m_m_i_n_g)', 'abs(X_B_l_a_c_k_m_a_n)', 'abs(X_c_h_e_b_w_i_n_1_0_0)', 'abs(X_c_h_e_b_w_i_n_1_2_0)');
xlabel('f = 0:0.1:500 Hz');



%----------Rysuje sygnaly i okna na jednym wykresie----------
figure(3);
subplot(3,1,1);
hold all;
plot(abs(x_r), 'r-');
plot(abs(x_h), 'g-');
plot(abs(x_b), 'c-');
plot(abs(x_c100), 'b-');
plot(abs(x_c120), 'k-');
title('Sygnaly po wymnozeniu przez okna');
legend('x_r_e_c_t_p_u_l_s', 'x_H_a_m_m_i_n_g', 'x_B_l_a_c_k_m_a_n', 'x_c_h_e_b_w_i_n_1_0_0', 'x_c_h_e_b_w_i_n_1_2_0');
xlabel('f = 0:0.1:500 Hz');
subplot(3,1,2);
hold all;
plot(rectpuls(t, (N-1)/f_s), 'r-');  %okno prostokatne
plot(hamming(N)', 'g-');             %okno Hamminga 
plot(blackman(N)', 'c-');            %okno Blackman
plot(chebwin(N, 100)', 'b-');        %okno Czebyszewa (100 dB)
plot(chebwin(N, 120)', 'k-');        %okno Czebyszewa (120 dB)
title('Okna przez ktore przemnozono oryginalny sygnal');
legend('okno prostokatne', 'okno Hamminga', 'okno Blackman', 'okno Czebyszewa (100 dB)', 'okno Czebyszewa (120 dB)');
xlabel('f = 0:0.1:500 Hz');
subplot(3,1,3)
plot(x, 'r-o');
title('Sygnal oryginalny x');
legend('x');
xlabel('f = 0:0.1:500 Hz');