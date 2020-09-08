function [E,freq,v] = PhaseShiftOfSW(rec,dt,offset,dx,vmin,dv,vmax,fmin,fmax)
[nt,nx] = size(rec);      % ��õ����¼��ά��
rec_fft = zeros(nt,nx);   % �����¼��Ƶ��

v=vmax:-dv:vmin;          % ���ٶ�ɨ������
lv = length(v);           % ���ٶ������ĳ���

Fs = 1/dt;                % ������
f = Fs*(0:(nt/2))/nt;     % Ƶ������

% ȷ��Ƶɢ����ͼ����СƵ����Ƶ�������е�����
ind = find(f>fmin);       
fmin_ind  = ind(1)-1;

% ȷ��Ƶɢ����ͼ�����Ƶ����Ƶ�������е�����
ind = find(f>fmax);
fmax_ind = ind(1);

freq = f(fmin_ind:fmax_ind); % Ƶɢ����ͼ�е�Ƶ������
lf = length(freq);

for j=1:nx
    rec_fft(:,j) = fft(rec(:,j),nt,1); % ��õ����¼��Ƶ��
end

rec_fft_amp = abs(rec_fft);

rec_fft_n = rec_fft./rec_fft_amp;      % �Ե����¼��Ƶ������һ��

E = zeros(lv,lf); % Ƶɢ��������


for j=1:lf
    f_ind = fmin_ind+j-1;
    for i=1:lv
        for k=1:nx
            x=offset+(k-1)*dx;
            E(i,j) = E(i,j)+exp(1i*2*pi*f(f_ind)*x/v(i))*rec_fft_n(f_ind,k);
        end
    end
    E(:,j) = abs(E(:,j)./max(abs(E(:,j))));  % ��ÿ��Ƶ�ʵ�Ƶɢ��������һ��
end
end