function [E,freq,v] = PhaseShiftOfSW(rec,dt,offset,dx,vmin,dv,vmax,fmin,fmax)
[nt,nx] = size(rec);      % 获得地震记录的维度
rec_fft = zeros(nt,nx);   % 地震记录的频谱

v=vmax:-dv:vmin;          % 相速度扫描向量
lv = length(v);           % 相速度向量的长度

Fs = 1/dt;                % 采样率
f = Fs*(0:(nt/2))/nt;     % 频率向量

% 确定频散能量图中最小频率在频率向量中的索引
ind = find(f>fmin);       
fmin_ind  = ind(1)-1;

% 确定频散能量图中最大频率在频率向量中的索引
ind = find(f>fmax);
fmax_ind = ind(1);

freq = f(fmin_ind:fmax_ind); % 频散能量图中的频率向量
lf = length(freq);

for j=1:nx
    rec_fft(:,j) = fft(rec(:,j),nt,1); % 获得地震记录的频谱
end

rec_fft_amp = abs(rec_fft);

rec_fft_n = rec_fft./rec_fft_amp;      % 对地震记录的频谱作归一化

E = zeros(lv,lf); % 频散能量矩阵


for j=1:lf
    f_ind = fmin_ind+j-1;
    for i=1:lv
        for k=1:nx
            x=offset+(k-1)*dx;
            E(i,j) = E(i,j)+exp(1i*2*pi*f(f_ind)*x/v(i))*rec_fft_n(f_ind,k);
        end
    end
    E(:,j) = abs(E(:,j)./max(abs(E(:,j))));  % 对每个频率的频散能量作归一化
end
end