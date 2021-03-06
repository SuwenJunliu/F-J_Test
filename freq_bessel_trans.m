function [I,f] = freq_bessel_trans(C,r,t,v)
%% transform (x-t) data into F-J data

%input: C       the cross-correlation result of station pair
%               2D matrix, size: (nt,nx)
%       r       the vector of offest
%               length: nx
%       t       the vector of time series, length:nt
%       v       the vector of velocity, length:nv
%       r0      the vector of C_r0 , length:nt

%output I       the omega-k result matrix. same size of C

% written by junliu  @ IGGCAS

% try to fix C_r0 ?


%% 
nx = length(r);
nt = length(t);

%new_r = zeros(1,nx+1);
%new_r(2:end) = r;
%r = new_r;

%new_C = zeros(nt,nx+1);
%new_C(:,1) = r0;
%new_C(:,2:end) = C ;
%C = new_C;


nv = length(v);
I = zeros(nv,nt);
pi = 3.1415927;
%fft_C_amp = abs(fft(C));
%fft_C = fft(C)./fft_C_amp;
fft_C = fft(C);
% complex number?
dt = abs(t(2) - t(1));

inum = 0:1:nt-1;
omega = pi/(dt*nt/2) * inum+0.00001;
% the first row is zero, so a/k equals INF?
% quick fix: plus 0.00001

f = inum /dt/nt;






% In fact, there is still a problem about r_0
% According to the doi:10.1029/2018JB016595 
% b_j = [C (r_j) - C(r_(j-1)) ]�M ��r_j
% when j = 1, what on earth C_r0 is ?
% For cross-correlation data, it might be the auto-correlation ?
% but for real data ?


% base on  doi:10.1029/2018JB016595 
for iOmega = 1:nt/2
    for iv = 1:nv
        temp = 0;
        for ir = 2:nx
            k = omega(iOmega) / v(iv);
            b = ((fft_C(iOmega,ir)-fft_C(iOmega,ir-1))/(r(ir)-r(ir-1)));
            a = fft_C(iOmega,ir-1) - r(ir-1) * b;
            
            %temp1 = a/k* ( r(ir)*besselj(0,k*r(ir)) - r(ir-1)*besselj(0,k*r(ir-1)) );
            %temp2 = b/k* ( r(ir).^2*besselj(1,k*r(ir)) - r(ir-1).^2*besselj(1,k*r(ir-1)) );
            %temp3 = b/k^2* ( (r(ir)*besselj(0,k*r(ir))) -r(ir-1)*besselj(0,k*r(ir-1)) );
            
            fun = @(bessel) besselj(0,k*bessel);
            %temp4 = - b/k^3* integral(fun,k*r(ir-1),k*r(ir));
            
            temp1 = a/k * ( fft_C(iOmega,ir) * r(ir) * besselj(1,k*r(ir)) - fft_C(iOmega,ir-1) * r(ir-1) * besselj(1,k*r(ir-1)));
            
            temp11 = b/k.^3 * ( (k*r(ir)).^2*besselj(1,k*r(ir)) - (k*r(ir-1)).^2*besselj(1,k*r(ir-1))  );
            temp2 = b/k.^3 * ( k*r(ir)*besselj(0,k*r(ir)) -  k*r(ir-1)*besselj(0,k*r(ir-1)));
            temp3 = b/k.^3 * integral(fun,r(ir-1),r(ir));
            
            temp = temp1+temp11 + temp2 - temp3;
           
            
                
            
        end
        I(iv,iOmega) = temp;
        
    end
    I(:,iOmega) = abs(I(:,iOmega)./max(abs(I(:,iOmega))));
    disp(['omega ' , num2str((iOmega/nt)*100),'%'])
end
   

end

