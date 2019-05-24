% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% Æwiczenie: Kompresja sygna³u mowy wed³ug standardu LPC-10
% ----------------------------------------------------------


%% Deklaracja danych
clear all; clf; close all; clc;

%[x,fpr,Nbits]=wavread('mowa1.wav');	      % wczytaj sygna³ mowy (ca³y)
[x,fpr]=audioread('mowa1.wav');     	      % wczytaj sygna³ mowy (ca³y)
%figure(1); plot(x); title('sygn³ mowy');      % poka¿ go
%soundsc(x,fpr);                               % oraz odtwó¿ go na g³oœnikach (s³uchawkach)

N=length(x);	  % d³ugoœæ sygna³u
Mlen=240;		  % d³ugoœæ okna Hamminga (liczba próbek)
Mstep=180;		  % przesuniêcie okna w czasie (liczba próbek)
Np=10;			  % rz¹d filtra predykcji
gdzie=Mstep+1;	  % pocz¹tkowe polo¿enie pobudzenia dŸwiêcznego

lpc=[];								% tablica na wspóczynniki modelu sygna³u mowy
s=[];                               % ca³a mowa zsyntezowana
ss=[];							    % fragment sygna³u mowy zsyntezowany
bs=zeros(1,Np);					    % bufor na fragment sygna³u mowy
Nramek=floor((N-Mlen)/Mstep+1);     % ile fragmentów (ramek) jest do przetworzenia


x=filter([1 -0.9735], 1, x);      % filtracja wstêpna (preemfaza) - opcjonalna

%% c) Progowanie sygna³u x:
    x_prog = zeros(length(x), 1);
    
    for n=1:length(x)
        if mod(n-1, Mstep) == 0 && n+Mlen<=length(x)         % Jeœli przeszliœmy do kolejnej ramki, obliczamy nowe P
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
    
    
%% Znajduje sygna³ resztkowy dla g³oski dŸwiêcznej "a"
for  nr = 1 : 19
    % pobierz kolejny fragment sygna³u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);                                   % usuñ wartoœæ œredni¹
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    
    offset=20; rmax=max( r(offset : Mlen) );        % znajdŸ maksimum funkcji autokorelacji
    imax=find(r==rmax);                             % znajdŸ indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end;  % g³oska dŸwiêczna/bezdŸwiêczna?
    
    if (T>80) T=round(T/2); end;                    % znaleziono drug¹ podharmoniczn¹
    rr(1:Np,1)=(r(2:Np+1))';
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];		% zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;									% oblicz wspóczynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;							% oblicz wzmocnienie
    H=freqz(1,[1;a]);								% oblicz jego odp. czêstotliwoœciow¹
    
    
    przedzial = ((nr-1)*Mstep+1):(nr*Mstep);                  % #1
    sygnal_resztkowy = filter(a, 1, x(przedzial));            % #2
    
end

T_res = T;
przedzial_res = przedzial;

disp(['Okres sygna³u resztkowego = ', num2str(T_res)]);

figure(1);
subplot(1,2,1);
plot(x(przedzial_res), 'b');
title('G³oska dŸwiêczna "a"');
xlabel('nr próbki');
subplot(1,2,2);
plot(sygnal_resztkowy, 'r');
title('Sygna³ resztkowy "a"');
xlabel('nr próbki');

figure(2);
plot(abs(H), 'k.-');
title('Cha-ka a-cz transmitancji flitru predykcji');
xlabel('[Hz]');









%% Analiza kolejnych ramek
for  nr = 1 : Nramek
    % pobierz kolejny fragment sygna³u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);                                   % usuñ wartoœæ œredni¹
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    
    offset=20; rmax=max( r(offset : Mlen) );        % znajdŸ maksimum funkcji autokorelacji
    imax=find(r==rmax);                             % znajdŸ indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end;  % g³oska dŸwiêczna/bezdŸwiêczna?
    
    if (T>80) T=round(T/2); end;                    % znaleziono drug¹ podharmoniczn¹
    rr(1:Np,1)=(r(2:Np+1))';
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];		% zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;									% oblicz wspóczynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;							% oblicz wzmocnienie
    H=freqz(1,[1;a]);								% oblicz jego odp. czêstotliwoœciow¹
    
    lpc=[lpc; T; wzm; a; ];						% zapamiêtaj wartoœci parametrów
    
    %% SYNTEZA - odtwórz na podstawie parametrów ----------------------------------------------------------------------
    if (T~=0) gdzie=gdzie-Mstep; end				% 'przenieœ' pobudzenie dŸwiêczne
    for n=1:Mstep
            if( T==0)
                pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1;			% pobudzenie szumowe
            else
                if (n==gdzie) pob=1; 
                    gdzie=gdzie+T;                 % pobudzenie dŸwiêczne
                else pob=0; end
            end
        
        %ss(n)=wzm*pob-bs*a;		% filtracja 'syntetycznego' pobudzenia
        ss(n)=wzm*sygnal_resztkowy(mod(n, T_res) + T_res)-bs*a; % gloska dziwieczna            % #3
        
        
                                % wzm*pob = err(n) = b³¹d predykcji (19.12b)
                                % -ba*a =  prognozowana wartoœæ sygna³u w chwili n-tej [PREDYKCJA] (19.12a)
        bs=[ss(n) bs(1:Np-1) ];	% przesuniêcie bufora wyjœciowego
    end
    s = [s ss];						% zapamiêtanie zsyntezowanego fragmentu mowy
    
end

s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

%% plot(s); title('mowa zsyntezowana'); pause
soundsc(s, fpr)
