[E,freq,v] = PhaseShiftOfSW(Uz(:,2:end),(t(2)-t(1)),100000,2500,3500,10,5500,0.1,7);

x=[freq(1) freq(end)];
y=[v(1) v(end)];

figure(1)
imagesc(x,y,E);
colormap(jet);
colorbar;
set(gca,'YDir','normal');
xlabel('频率/(Hz)');
ylabel('相速度/(m/s)');