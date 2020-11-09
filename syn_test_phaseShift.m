% SEIS2
dt = 0.001;
offset = 0.5;
dx = 0.002;
[E,freq,v] = PhaseShiftOfSW(suma_sg2,dt,0.001,dx,0.5,0.005,1.5,10,50);
imagesc(freq,v,E);
hold on
plot(data(:,1),data(:,3))
set(gca,'YDir','normal')

%SEIS4
%dt = 0.0195;
%offset = 0.500;
%dx = 0.505;
%[E,freq,v] = PhaseShiftOfSW(seis4(27:end,30:end),dt,offset,dx,2.5,0.005,5.5,0.01,10);



[Iq,fq] = freq_bessel_trans(suma_sg2,x(1:50)/1000,t_total/1000,flip(v));