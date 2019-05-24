%Lab05 -  Implementacja algorytmu radix-2
%mag2db(y) = 20log(y);

close all;
clear all;
clc;

%============================Dane==============================
N = 1024;           %Dlugosc sygnalu
f_s = 1000;         %czestotliwosc probkowania [Hz]

%==========Generowanie wektora zespolonego o rozkladzie normalnym==========
x_real = zeros(1, N);
x_imaginary = zeros(1, N);
for i=1:N
    x_real(i) = normrnd(0, 1);   %Czesc rzeczywista sygnalu losowego o rozkladzie normalnym
    x_imaginary(i) = normrnd(0, 1);   %Czesc urojona sygnalu losowego o rozkladzie normalnym
end


%===================Zapis wektora do plikow====================
x = [x_real; x_imaginary]';
save('x.mat', 'x', '-ascii');
save('xcpp.dat','x','-ascii','-double','-tabs');








%   double tab[1024][2];
%	ifstream xcc_dat;
%	xcc_dat.open("C:\\Users\\pawel\\Documents\\MATLAB\\CPS_Lab04_solved\\xcpp.dat");
%	for (int i = 0; i < 1024; ++i) {
%		for (int j = 0; j < 2; ++j) {
%			xcc_dat >> tab[i][j];
%		}
%	}