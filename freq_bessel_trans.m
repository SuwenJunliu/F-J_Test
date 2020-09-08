function [I,f] = freq_bessel_trans(C,r,t,v)
%% transform (x-t) data into F-J data

%input: C       the cross-correlation result of station pair
%               2D matrix, size: (nt,nx)
%       r       the vector of offest
%               length: nx
%       t       the vector of time series, length:nt
%       v       the vector of velocity, length:nv

%output I       the omega-k result matrix. same size of C

% written by junliu 2020/9/4 @ IGGCAS


%% 

nx = length(r);
nt = length(t);
nv = length(v);
I = zeros(nt,nv);
pi = 3.1415927;
fft_C = fft(C);
% complex number?
dt = abs(t(2) - t(1));

inum = 0:1:nt-1;
omega = pi/(dt*nt/2) * inum+0.00001;
% the first row is zero, so a/k equals INF?
% quick fix: plus 0.00001

f = inum /dt/nt;



% this part base on "背景噪音提取高阶频散曲线的矢量波数变换方法" cnki

for iOmega = 1:nt
    for iv = 1:nv
        temp = 0;
        for ir = 2:nx
            k = omega(iOmega) / v(iv);
            b = ((fft_C(iOmega,ir)-fft_C(iOmega,ir-1))/(r(ir)-r(ir-1)));
            a = fft_C(iOmega,ir-1) - r(ir-1) * b;
            
            temp1 = a/k* ( r(ir)*besselj(1,k*r(ir)) - r(ir-1)*besselj(1,k*r(ir-1)) );
            temp2 = b/k^2* ( k*r(ir)^2*besselj(1,k*r(ir)) - k*r(ir-1)^2*besselj(1,k*r(ir-1)) );
            temp3 = b/k^2* ( (r(ir)*besselj(0,k*r(ir))) -r(ir-1)*besselj(0,k*r(ir-1)) );
            
            fun = @(bessel) besselj(0,k*bessel);
            temp4 = - b/k^2* integral(fun,r(ir-1),r(ir));
            
            % sometimes negative?
            temp = temp1 + temp2 + temp3 + temp4;
                
            
        end
        I(iOmega,iv) = temp;
        
    end
    disp(['omega ' , num2str((iOmega/nt)*100),'%'])
end
   

end

