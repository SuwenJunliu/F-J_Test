%% Frequency-Bessel Method 

r0 = 50;
dz0 = 2.5;
inum = 50;
x = r0+dz0:dz0:r0+dz0*inum;
x = x * 1000;
v = linspace(3.5,5.5,100);
v = v * 1000;

tmax = 140;
t = Uz(:,1);

data = Uz(:,2:end);
r0 = zeros(1,2048);

[I,f] = freq_bessel_trans(data,x,t,v,r0);