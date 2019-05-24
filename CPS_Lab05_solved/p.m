function wynik = p(N)
    w_3dB = 2*pi*100;   %[rd/s]
    k = 1:N;
    wynik = w_3dB * exp( 1i*( pi/2 + pi/( 2*N) + (k-1)*pi/N));
end