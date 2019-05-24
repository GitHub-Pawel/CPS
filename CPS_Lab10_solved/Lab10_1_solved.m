% ----------------------------------------------------------
% Tabela 19-4 (str. 567)
% �wiczenie: Kompresja sygna�u mowy wed�ug standardu LPC-10
% ----------------------------------------------------------


%% Deklaracja danych
clear all; clf; close all; clc;

%[x,fpr,Nbits]=wavread('mowa1.wav');	      % wczytaj sygna� mowy (ca�y)
[x,fpr]=audioread('mowa1.wav');     	      % wczytaj sygna� mowy (ca�y)
figure(1); plot(x); title('sygn� mowy');      % poka� go
soundsc(x,fpr);                               % oraz odtw� go na g�o�nikach (s�uchawkach)

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

%% 1. Wyb�r fragmet�w mowy, zawieraj�cych odpowiednie g�oski:

gloska_dzw = 3200:3600;    % g�oska "a" - mowa 1
%gloska_dzw = 183000:183240;    % g�oska d�wi�czna - mowa 2
%gloska_dzw = 11700:11940;    % g�oska d�wi�czna - mowa 3
%soundsc(x(gloska_dzw),fpr);            % ODTW�RZ

gloska_bez = 37460:37860;    % g�oska "cz" - mowa 1
%gloska_bez = 6300:6540;    % g�oska bezd�wi�czna - mowa 2
%gloska_bez = 12060:12300;    % g�oska bezd�wi�czna - mowa 2
%soundsc(x(gloska_bez),fpr);            % ODTW�RZ

stan_przej = 3670:4070;      %stan przej�ciowy "a" -> "t"   - mowa 1
%stan_przej = 4500:4740;      %stan przej�ciowy - mowa 2
%stan_przej = 11880:12120;      %stan przej�ciowy - mowa 3
%soundsc(x(stan_przej), fpr);           % ODTW�RZ

%% a) sygna� czasowy oraz widmo g�sto�ci widmowej mocy sygna�u przed oraz po preemfazie,

figure(2);
subplot(3,1,1);
plot(x(gloska_dzw), 'k.-');
title('g�oska d�wi�czna "a"');
xlabel("nr pr�bki");
subplot(3,1,2);
plot(x(gloska_bez), 'k.-');
title('g�oska bezd�wi�czna "cz"');
xlabel("nr pr�bki");
subplot(3,1,3);
plot(x(stan_przej), 'k.-');
title('stan przej�ciowy "a...t"');
xlabel("nr pr�bki");


figure(3);
subplot(3,1,1);
plot(dspdata.psd(abs(x(gloska_dzw)),'Fs',fpr));
title('Power Spectral Density: "a" [przed preemfaz�]');
subplot(3,1,2);
plot(dspdata.psd(abs(x(gloska_bez)),'Fs',fpr));
title('Power Spectral Density: "cz" [przed preemfaz�]');
subplot(3,1,3);
plot(dspdata.psd(abs(x(stan_przej)),'Fs',fpr));
title('Power Spectral Density: "a...t" [przed preemfaz�]');


x=filter([1 -0.9735], 1, x);      % filtracja wst�pna (preemfaza) - opcjonalna

figure(2);
subplot(3,1,1);
hold all;
plot(x(gloska_dzw), 'c.-');
legend('[przed preemfaz�]', '[po preemfazie]');
xlabel("nr pr�bki");
subplot(3,1,2);
hold all;
plot(x(gloska_bez), 'c.-');
legend('[przed preemfaz�]', '[po preemfazie]');
xlabel("nr pr�bki");
subplot(3,1,3);
hold all;
plot(x(stan_przej), 'c.-');
legend('[przed preemfaz�]', '[po preemfazie]');
xlabel("nr pr�bki");


figure(4);
subplot(3,1,1);
plot(dspdata.psd(abs(x(gloska_dzw)),'Fs',fpr));
title('Power Spectral Density: "a" [po preemfazie]');
subplot(3,1,2);
plot(dspdata.psd(abs(x(gloska_bez)),'Fs',fpr));
title('Power Spectral Density: "cz" [po preemfazie]');
subplot(3,1,3);
plot(dspdata.psd(abs(x(stan_przej)),'Fs',fpr));
title('Power Spectral Density: "a...t" [po preemfazie]');

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
    
    figure(7);
    subplot(3,1,1);
    hold all;
    plot(x(gloska_dzw), 'b.-');
    plot(x_prog(gloska_dzw), 'r.-');
    title('Gloska dzwieczna "a"');
    legend('Przed progowaniem', 'Po progowaniu');
    xlabel('Nr probki');
    subplot(3,1,2);
    hold all;
    plot(x(gloska_bez), 'b.-');
    plot(x_prog(gloska_bez), 'r.-');
    title('Gloska bezdzwieczna "cz"');
    legend('Przed progowaniem', 'Po progowaniu');
    xlabel('Nr probki');
    subplot(3,1,3);
    hold all;
    plot(x(stan_przej), 'b.-');
    plot(x_prog(stan_przej), 'r.-');
    title('Stan przej�ciowy "a...t"');
    legend('Przed progowaniem', 'Po progowaniu');
    xlabel('Nr probki');
    
    %x(:) = x_prog(:);
%%


%% Analiza kolejnych ramek

for  nr = 1 : Nramek
    clc;
    disp(['Ramka nr: ', num2str(nr), ' zakres: ', num2str((nr-1)*Mstep), '-', num2str((nr-1)*Mstep+Mlen)]);
    soundsc(x(((nr-1)*Mstep+1):((nr-1)*Mstep+Mlen)+1),fpr);
    
    % pobierz kolejny fragment sygna�u
    n = 1+(nr-1)*Mstep : Mlen + (nr-1)*Mstep;
    bx = x(n);
    
    % ANALIZA - wyznacz parametry modelu ---------------------------------------------------
    bx = bx - mean(bx);                                   % usu� warto�� �redni�
    for k = 0 : Mlen-1
        r(k+1) = sum( bx(1 : Mlen-k) .* bx(1+k : Mlen) ); % funkcja autokorelacji
    end
    figure(8); subplot(411); plot(n,bx); title('fragment sygna�u mowy');
    figure(8); subplot(412); plot(r); title('jego funkcja autokorelacji');
    hold on; 
    [ymax2,imax2,ymin2,imin2] = extrema(r);
    for i=1:length(ymax2)
        if ymax2(i) > 0.35*r(1)
            plot(imax2(i),ymax2(i),'r*');
        end
    end
    hold off;
    
    offset=20; rmax=max( r(offset : Mlen) );        % znajd� maksimum funkcji autokorelacji
    imax=find(r==rmax);                             % znajd� indeks tego maksimum
    if ( rmax > 0.35*r(1) ) T=imax; else T=0; end;  % g�oska d�wi�czna/bezd�wi�czna?
    if (T>80) T=round(T/2); end;                    % znaleziono drug� podharmoniczn�
    disp(['Okres T= ', num2str(T)]);				% wy�wietl warto�� T
    rr(1:Np,1)=(r(2:Np+1))';
    for m=1:Np
        R(m,1:Np)=[r(m:-1:2) r(1:Np-(m-1))];		% zbuduj macierz autokorelacji
    end
    a=-inv(R)*rr;									% oblicz wsp�czynniki filtra predykcji
    wzm=r(1)+r(2:Np+1)*a;							% oblicz wzmocnienie
    H=freqz(1,[1;a]);								% oblicz jego odp. cz�stotliwo�ciow�
    %%  b) charakterystyk� amplitudowo-cz�stotliwo�ciow� filtru H(z):
    
    if nr == 19
        figure(5);
        subplot(3,1,1);
        plot(abs(H), 'b.-');
    elseif nr == 20
        figure(5);
        hold all;
        plot(abs(H), 'r-.');
        title('Cha-ka a-cz transmitancji dla "a"');
        legend('3240-3360', '3360-3600');
        xlabel('[Hz]');
    end
    
    
    if nr == 209
        figure(5);
        subplot(3,1,2);
        plot(abs(H), 'b.-');
    elseif nr == 210
        figure(5);
        hold all;
        plot(abs(H), 'r-.');
        title('Cha-ka a-cz transmitancji dla "cz"');
        legend('37440-37620', '37620-17800');
        xlabel('[Hz]');
    end
    
    
    if nr == 22
        figure(5);
        subplot(3,1,3);
        plot(abs(H), 'b.-');
    elseif nr == 23
        figure(5);
        hold all;
        plot(abs(H), 'r-.');
        title('Cha-ka a-cz transmitancji dla "a...t"');
        legend('3600-3780', '3780-3960');
        xlabel('[Hz]');
    end
   
    
    figure(8); subplot(413); plot(abs(H)); title('widmo filtra traktu g�osowego');
    hold on; 
    [ymax,imax,ymin,imin] = extrema(abs(H));
    plot(imax,ymax,'r*');
    hold off;
    disp('Czestotliwo�ci tonu podstawowowego: ');
    for i=1:length(imax)
        disp([num2str(i), '. ', num2str(imax(i)), ' [Hz]']);
    end
    
    lpc=[lpc; T; wzm; a; ];						% zapami�taj warto�ci parametr�w
    
    % SYNTEZA - odtw�rz na podstawie parametr�w ----------------------------------------------------------------------
    % T = 0;                                        % usu� pierwszy znak '%' i ustaw: T = 80, 50, 30, 0 (w celach testowych)
    if (T~=0) gdzie=gdzie-Mstep; end				% 'przenie�' pobudzenie d�wi�czne
    for n=1:Mstep
        % T = 70; % 0 lub > 25 - w celach testowych
        if( T==0)
            pob=2*(rand(1,1)-0.5); gdzie=(3/2)*Mstep+1;			% pobudzenie szumowe
            pobudzenie = 'pobudzenie: szumowe';
        else
            if (n==gdzie) pob=1; 
                gdzie=gdzie+T;                 % pobudzenie d�wi�czne
                pobudzenie = 'pobudzenie: d�wi�czne';
            else pob=0; end
        end
        ss(n)=wzm*pob-bs*a;		% filtracja 'syntetycznego' pobudzenia
        bs=[ss(n) bs(1:Np-1) ];	% przesuni�cie bufora wyj�ciowego
    end
    disp(pobudzenie);
    figure(8); subplot(414); plot(ss); title('zsyntezowany fragment sygna�u mowy'); %pause
    s = [s ss];						% zapami�tanie zsyntezowanego fragmentu mowy
    
end

s=filter(1,[1 -0.9735],s); % filtracja (deemfaza) - filtr odwrotny - opcjonalny

%% d) Por�wnaj ramk� oryginaln� i zsyntezowan� w dziedzinie czasu oraz cz�stotliwo�ci:

figure(9);
subplot(3,1,1);
plot(x(gloska_dzw), 'k.-');
hold all;
plot(s(gloska_dzw), 'r.-');
title('g�oska d�wi�czna "a"');
xlabel("nr pr�bki");
legend('oryginalny', 'zsyntezowany');
subplot(3,1,2);
plot(x(gloska_bez), 'k.-');
hold all;
plot(s(gloska_bez), 'r.-');
title('g�oska bezd�wi�czna "cz"');
xlabel("nr pr�bki");
legend('oryginalny', 'zsyntezowany');
subplot(3,1,3);
plot(x(stan_przej), 'k.-');
hold all;
plot(s(stan_przej), 'r.-');
title('stan przej�ciowy "a...t"');
xlabel("nr pr�bki");
legend('oryginalny', 'zsyntezowany');


figure(10);
subplot(3,1,1);
plot(abs(fft(x(gloska_dzw))), 'k.-');
hold all;
plot(abs(fft(s(gloska_dzw))), 'r.-');
title('Cha-ka a-cz: g�oska d�wi�czna "a"');
xlabel("[Hz]");
legend('oryginalny', 'zsyntezowany');
subplot(3,1,2);
plot(abs(fft(x(gloska_bez))), 'k.-');
hold all;
plot(abs(fft(s(gloska_bez))), 'r.-');
title('Cha-ka a-cz: g�oska bezd�wi�czna "cz"');
xlabel("[Hz]");
legend('oryginalny', 'zsyntezowany');
subplot(3,1,3);
plot(abs(fft(x(stan_przej))), 'k.-');
hold all;
plot(abs(fft(s(stan_przej))), 'r.-');
title('Cha-ka a-cz: stan przej�ciowy "a...t"');
xlabel("[Hz]");
legend('oryginalny', 'zsyntezowany');

%% plot(s); title('mowa zsyntezowana'); pause
soundsc(s, fpr)
