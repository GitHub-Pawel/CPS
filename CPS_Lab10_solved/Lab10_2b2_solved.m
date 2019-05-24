% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% �wiczenie: Kompresja sygna�u mowy wed�ug standardu LPC-10
% ----------------------------------------------------------


%% Deklaracja danych
clear all; clf; close all; clc;

%[x,fpr,Nbits]=wavread('mowa1.wav');	      % wczytaj sygna� mowy (ca�y)
[x,fpr]=audioread('mowa1.wav');     	      % wczytaj sygna� mowy (ca�y)
%figure(1); plot(x); title('sygn� mowy');      % poka� go
%soundsc(x,fpr);                               % oraz odtw� go na g�o�nikach (s�uchawkach)

N=length(x);	  % d�ugo�� sygna�u
Mlen=240;		  % d�ugo�� okna Hamminga (liczba pr�bek)
Mstep=180;		  % przesuni�cie okna w czasie (liczba pr�bek)
Np=10;			  % rz�d filtra predykcji
gdzie=Mstep+1;	  % pocz�tkowe polo�enie pobudzenia d�wi�cznego

lpc=[];								% tablica na wsp�czynniki modelu sygna�u mowy
s=[];                               % ca�a mowa zsyntezowana
ss=[];							    % fragment sygna�u mowy zsyntezowany
bs=zeros(1,Np);					    % bufor na fragment sygna�u mowy
Nramek=floor((N-Mlen)/Mstep+1);     % ile fragment�w (ramek) jest do przetworzenia


x=filter([1 -0.9735], 1, x);      % filtracja wst�pna (preemfaza) - opcjonalna

%% c) Progowanie sygna�u x:
    x_prog = zeros(length(x), 1);
    
    for n=1:length(x)
        if mod(n-1, Mstep) == 0 && n+Mlen<=length(x)         % Je�li przeszli�my do kolejnej ramki, obliczamy nowe P
            P = 0.3*max(x(n:n+Mlen));                        % Progujemy kazda ramke z osobna?
        end
        
        if x(n) >= P
            x_prog(n) = x(n) - P;
        elseif x(n) <= -P
            x_prog(n) = x(n) + P;
        else
            x_prog(n) = 0;
        end
    end
    
    %x(:) = x_prog(:);
    
    
%% Znajduje sygna� resztkowy dla g�oski d�wi�cznej "a"
for  nr = 1 : 19
    % pobierz kolejny fragment sygna�u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);                                   % usu� warto�� �redni�
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    
    offset=20; rmax=max( r(offset : Mlen) );        % znajd� maksimum funkcji autokorelacji
    imax=find(r==rmax);                             % znajd� indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end;  % g�oska d�wi�czna/bezd�wi�czna?
    
    if (T>80) T=round(T/2); end;                    % znaleziono drug� podharmoniczn�
    rr(1:Np,1)=(r(2:Np+1))';
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];		% zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;									% oblicz wsp�czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;							% oblicz wzmocnienie
    H=freqz(1,[1;a]);								% oblicz jego odp. cz�stotliwo�ciow�
    
    
    przedzial = ((nr-1)*Mstep+1):(nr*Mstep);                  % #1
    sygnal_resztkowy = filter(a, 1, x(przedzial));            % #2
    
end

T_res = T;
przedzial_res = przedzial;

disp(['Okres sygna�u resztkowego = ', num2str(T_res)]);

figure(1);
subplot(1,2,1);
plot(x(przedzial_res), 'b');
title('G�oska d�wi�czna "a"');
xlabel('nr pr�bki');
subplot(1,2,2);
plot(sygnal_resztkowy, 'r');
title('Sygna� resztkowy "a"');
xlabel('nr pr�bki');

figure(2);
plot(abs(H), 'k.-');
title('Cha-ka a-cz transmitancji flitru predykcji');
xlabel('[Hz]');









%% Analiza kolejnych ramek
for  nr = 1 : Nramek
    % pobierz kolejny fragment sygna�u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);                                   % usu� warto�� �redni�
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    
    offset=20; rmax=max( r(offset : Mlen) );        % znajd� maksimum funkcji autokorelacji
    imax=find(r==rmax);                             % znajd� indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end;  % g�oska d�wi�czna/bezd�wi�czna?
    
    if (T>80) T=round(T/2); end;                    % znaleziono drug� podharmoniczn�
    rr(1:Np,1)=(r(2:Np+1))';
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];		% zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;									% oblicz wsp�czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;							% oblicz wzmocnienie
    H=freqz(1,[1;a]);								% oblicz jego odp. cz�stotliwo�ciow�
    
    lpc=[lpc; T; wzm; a; ];						% zapami�taj warto�ci parametr�w
    
    %% SYNTEZA - odtw�rz na podstawie parametr�w ----------------------------------------------------------------------
    if (T~=0) gdzie=gdzie-Mstep; end				% 'przenie�' pobudzenie d�wi�czne
    for n=1:Mstep
            if( T==0)
                pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1;			% pobudzenie szumowe
            else
                if (n==gdzie) pob=1; 
                    gdzie=gdzie+T;                 % pobudzenie d�wi�czne
                else pob=0; end
            end
        
        %ss(n)=wzm*pob-bs*a;		% filtracja 'syntetycznego' pobudzenia
        ss(n)=wzm*sygnal_resztkowy(mod(n, T_res) + T_res)-bs*a; % gloska dziwieczna            % #3
        
        
                                % wzm*pob = err(n) = b��d predykcji (19.12b)
                                % -ba*a =  prognozowana warto�� sygna�u w chwili n-tej [PREDYKCJA] (19.12a)
        bs=[ss(n) bs(1:Np-1) ];	% przesuni�cie bufora wyj�ciowego
    end
    s = [s ss];						% zapami�tanie zsyntezowanego fragmentu mowy
    
end

s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

%% plot(s); title('mowa zsyntezowana'); pause
soundsc(s, fpr)
