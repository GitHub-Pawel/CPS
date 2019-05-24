clear all;
clc;

load('adsl_x.mat');

plot(xcorr(x(:)));

j=1;
for i=1:2049
    i_backup = i;
    licznik = 1;
   if(abs(x(j)-x(i))<0.1)
       while(j<2049 && i<2049 && abs(x(j)-x(i))<0.1)
            j = j+1;
            i = i+1;
            licznik = licznik + 1;
       end
       if (licznik > 1)
            disp(licznik);
       end
   end
   j = 1;
   i = i_backup;
end