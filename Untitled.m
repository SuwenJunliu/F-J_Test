%% Frequency-Bessel Method 

r0 = 100;
dz0 = 5;
inum = 20;
x = r0+dz0:dz0:r0+dz0*inum;
x = x * 1000;
v = linspace(3.5,5.5,100);
v = v * 1000;

tmax = 140;
t = Urnew(:,1);

data = Urnew(:,2:end);

[I,f] = freq_bessel_trans(data(1:1000,:),x,t(1:1000),v);