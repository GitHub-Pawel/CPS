close all;
clear all;
clc;

N=20;
s(1:N)=sqrt(2/N);   %s_0
s(1)=sqrt(1/N);     %s_k, k~=0

A=zeros(N);     
for k=1:N           %kolejne wiersze macierzy
    for n=1:N       %kolejne kolumny macierzy
        A(k,n) = s(k) * cos(pi * (k-1)/N * (n-1 + 0.5));
                    %k oraz n nale¿y przesun¹æ o 1, aby k=n=0...N-1
    end
end

%Sprawdzam, czy wiersze macierzy A sa ortogonalne [metody z polecenia]
count = 0;
for i=1:N-1
    for j=i+1:N
        count = count+1;
        
        prod1 = sum(A(i,:) .* A(j,:));  % Iloczyn odpowiadaj¹cych sobie próbek
        prod2 = dot(A(i,:), A(j,:));    % w przestrzeni Euklidesowej
        prod3 = (A(i,:) * A(j,:)');     % bezpoœrednie obliczenie (mno¿enie wektorowe)
                                        % ,,0'' oznacza ¿e wektory s¹ ortogonalne 
        
        if(abs(prod1) > 1e-12 || abs(prod2) > 1e-12 || abs(prod3) > 1e-12)
            disp('Blad! Wiersze nie sa ortogonalne!');
        end
        
        %disp([num2str(count), '. ', num2str(prod1)]);
    end
end





%Sprawdzam ortogonalnosc/ortonormalnosc macierzy A [wzory 2.68/2.69]
ortogonalna = 1;
ortonormalna = 1;
for i=1:N
    for j=i:N
        if (i~=j &&  abs(sum(A(i,:).*conj(A(j,:)))) > 1e-12)
            ortogonalna = 0;
        elseif (i==j &&  abs(sum(A(i,:).*conj(A(j,:))) - 1) > 1e-12)
            ortonormalna = 0;
        end
    end
end

if (ortogonalna == 1 && ortonormalna == 1)
    disp('Macierz A jest ortonormalna');
elseif (ortogonalna == 1)
    disp('Macierz A jest ortogonalna');
else
    disp('Macierz A nie jest ortogonalna');
end
