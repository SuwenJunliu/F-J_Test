dt = 0.041;
offset = 0.5;
dx = 0.5;
[E,freq,v] = PhaseShiftOfSW(seis2,dt,offset,dx,2.5,0.005,5.5,0.01,10);
imagesc(E);