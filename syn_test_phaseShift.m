% SEIS2
dt = 0.041;
offset = 0.5;
dx = 0.5;
[E,freq,v] = PhaseShiftOfSW(seis2(27:end,2:end),dt,offset,dx,2.5,0.005,5.5,0.01,10);
imagesc(v,freq,E);

%SEIS4
dt = 0.0195;
offset = 0.500;
dx = 0.505;
[E,freq,v] = PhaseShiftOfSW(seis4(27:end,30:end),dt,offset,dx,2.5,0.005,5.5,0.01,10);