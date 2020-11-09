% DISECA ?  a Matlab code for dispersive waveform calculation, created by Renata Gazdova (gazdova@irsm.cas.cz) and Jan Vilhelm
% (vilhelm@natur.cuni.cz)

%Reference: Gaždov?, R., Vilhelm, J., 2011: DISECA ? A Matlab code for dispersive waveform calculations. 
%           Computers and Geotechnics, Vol. 38, No. 4, 526?531, doi:10.1016/j.compgeo.2011.03.001.

% Usage:
% DISECA: Program for model calculation of dispersive waveform. Calculation is based on summation the particular frequency 
% components with the introduction of phase shifts, corresponding to velocity dispersion and given distance from a seismic source.
% DISECA program computes the synthetic seismogram of a single mode dispersive wave without interference with other waves, and 
% with an exactly-defined dispersion curve, even over a very wide range of frequencies. 
% The knowledge of physical parameters of the medium is not required. The resulting waveform only contains an individual 
% dispersive wave of the selected mode, thus being particularly suitable for testing of methodologies for dispersive wave analysis.

%List of required input parameters:
% --dispersion curve of phase velocity (dependence of velocity on frequency) - two columns, with frequency in one and phase velosity in the second 
%                                       - the dispersion curve can represent an arbitrary mode of any dispersive wave, line 67,68
% geometry of the modeled synthetic wavefield: 
%        --maximal offset (maximal distance between the source and the receivers), line 56
%        --distance between receivers, line 57
% time parameters of the modeled synthetic waveform:
%          --duration of the seismic impulse F(t), line 58       
%          --sampling frequency - in time, line 59
%          --duration of the whole synthetic waveform, line 60

%     
% -- gain of the synthetic waveform, if it is not is required, set Z = 1, line 62
% -- damping coefficient - the damping coefficient ? represent decrease of the individual frequency component in amplitude with 
%                          time and is connected with absorption of seismic waves [23]. The damping coefficient ? can be chosen 
%                          as constant, linearly dependent on the time t, or with some more specific dependence on time. It can 
%                          also be chosen as frequency dependent (for example higher frequencies can have higher damping 
%                          coefficient than the lower ones). These choices can strongly affect the shape of the resulting 
%                          waveform, and can be used to adjust the fit of the measured data in greater detail. 
%             -- see line 62, 79-82

%List of optional input parameters:
% -- amplitude response - the amplitude represents the amplitude of the individual frequency component, line 69
% -- efect of geometrical spreading - it is possible to use by uncomment line 135

%List of output parameters:
% -- the final synthetic wavefiled is saved in matlab file: SyntheticWaveForm.m
%         -- file include: 'suma_sg' - Matrix, where each column represent values of amplitude of the single synthetic waveform 
%                          't_total' - vector that represent time
%                          'samplingInterval' - value of sampling interval
%                          'x - vector that represent distances (offsets} 

% Copyright 2008, Renata Gazdova and Jan Vilhlem; gazdova@irsm.cas.cz

% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details:
% http://www.gnu.org/licenses/.

clear all
close all

%parameters that need to be set, DISTANCES = METERS [m], TIME = MILISECONDS [ms]
max_x = 101;                        %maximal offset [m]
krok_x = 2;                         %distance between receivers [m] 
delka_impulsu = 300;                %duration of the seismic impulse F(t), in [ms]
samplingInterval = 1;               %sampling frequency - in time [ms]
delka_t = 1000;                     %duration of the whole synthetic waveform  [ms]
Z = 5;                              %gain of the synthetic waveform
beta = 1/8000;                      %damping coefficient

%load data of dispersive curve (in EXCEL format) - need to be defined in which columns are the frequency and the phase celocity data
data = xlsread('DispData_INPUT.xls');           %basic mode
%data = xlsread('DispData_INPUT_1mode.xls');    %first higher mode
f = data(:,1);                                  %column with frequency [Hz]
vph = data(:,3);                                %column with phase velocity [km/s]
%a = data(:,5);                                 %colums with AMPLITUDE RESPONSE - if we had this information
pocet = length(f);                              %calculation of the number of frequencz components
omega(1:pocet,1) = 2*pi*f*(10^-3);              %recalculation of frequency into the angular frequency 

x = 1:krok_x:max_x;                              %counting of offset axes x
max_k = max_x/krok_x;                            %counting of number of cycles that need to be done

%calculation of the seisic impulse for certain frequency components F(t)
t = 0:samplingInterval:delka_impulsu;              %time v ms    
    for i = 1:pocet      
        %betaa(i,1) = 1/(50*exp(1/f(i,1)*15));                  %exponencial increase of damping coefficient with frequency
        betaa(i,1) = (1*f(i,1))/8000;                           %linear increase of damping coefficient with frequency
        sg(i,:) = (exp(-1*beta*t)).*sin(t*omega(i,1));          %constat damping ratio
        sg2(i,:) = (exp(-1*betaa(i,1)*t)).*sin(t*omega(i,1));   %variable damping ratio
    end

    %plot of the all frequency components in case of variable damping ratio
    figure
    for i=1:pocet
        max_trace=max(sg(i,:));
        min_trace=min(sg(i,:));
        max_ampl_trace=max([max_trace, abs(min_trace)]);
        plot(t,sg2(i,:)*2+f(i,:));
        %plot(t,sg(i,:)/max_ampl_trace+f(i,:));
        hold on
    end
    hold off
    title('variable damping ratio')
    xlabel('time [ms]')
    ylabel('frequency component [Hz]')
    
    %plot of the all frequency components in case of constant damping ratio
    figure
    for i=1:pocet
        max_trace=max(sg(i,:));
        min_trace=min(sg(i,:));
        max_ampl_trace=max([max_trace, abs(min_trace)]);
        plot(t,sg(i,:)*2+f(i,:));
        %plot(t,sg(i,:)/max_ampl_trace+f(i,:));
        hold on
    end
    hold off
    title('constant damping ratio')
    xlabel('time [ms]')
    ylabel('frequency component [Hz]')

figure
%  1) calculation of the synthetic waveform from the dispersive curve - FIRST cycle for calculation F(t,f) for different offsets x
for k = 1:(max_x/krok_x)
    t_start(1:pocet,k) = x(1,k)./vph(:,1); %calculation of the arrival time of the given frequency component into the distance x [m]

    %calculation of the time offset - each frequency component have different velocity => different arrival time
    t_total = 0:samplingInterval:delka_t;
    vel_z = delka_t*(1/samplingInterval)+1;
    z = zeros(pocet,vel_z);
    t_start2(1:pocet,k) = round(t_start(:,k)*(1/samplingInterval)); 
    
    for i = 1:pocet    
        z_start(i,1) = t_start2(i,k);    
        z_end(i,1) = t_start2(i,k)+delka_impulsu*(1/samplingInterval); 
        %z(i,z_start(i,1):z_end(i,1))=sg(i,:);    %constant damping ratio
        z(i,z_start(i,1):z_end(i,1))=sg2(i,:);    %variable damping ratio
    end
    
    %sumation of the time shofted signals of different frequencies
    suma_sg(:,k) = sum(z);           %without geometrical spreading   
    %suma_sg(:,k) = 1/x(:,k)*sum(z);   %with geometrical spreading (the decrease of the energy with distance) 
    %plot(t_total,suma_sg(:,k)+k*krok_x)     %without normalization of the amplitude
    max_ampl=max(suma_sg(:,k));
    min_ampl=min(suma_sg(:,k));
    max_plot=max([max_ampl, abs(min_ampl)]);
    
    subplot(1,2,1)
    plot(t_total,suma_sg(:,k)*Z/max_plot+k*krok_x)   
    title('dispersive waveform')
    xlabel('time [ms]')
    ylabel('offset [m]')
    hold on
end
hold off
axis tight
%xlim([0 500])

%SPECTRUM
[vzorky,geofony]=size(suma_sg);
maxf = 1/(2*samplingInterval/1000);  
krokf=1/(t_total(end)/1000);         
frekv=0:krokf:maxf;
prumery_SumaSg = mean(suma_sg);
prumer1 = mean(prumery_SumaSg);
suma_sg2 = suma_sg-prumer1;
spe_z=abs(fft(suma_sg2));¡¤
spe_z2 = spe_z(1:((vzorky+1)/2),:);  

%normalization of the spectrum
for i=1:geofony
    n_spe_z2(:,i)=spe_z2(:,i)/max(spe_z(:,i));
end

%define frequency range of spectrum that will be drawn
fstart = 2;   %[Hz], need to be even number
fend = 180;   %[Hz], need to be even number
ZPS = 10 ;     %gain of the power spectrum

%figure
subplot(1,2,2)
for i=1:geofony
    plot(frekv,n_spe_z2(:,i)*ZPS+i);                                                         %plot all
    %plot(frekv,spe_z2(:,i)+i);                                                              %plot without normalization
    %plot(frekv(fstart/krokf:fend/krokf),n_spe_z2(fstart/krokf:fend/krokf,i*ZPS+i*krok_x);   %plot of the choosen part
    hold on
end
%title('amplitudove frekvencni spektrum, beta = linearne zavisla');
xlabel('frequency [Hz]');
%xlim([0 100]);
ylabel('offset [m]');
hold off
axis tight

%save('SyntheticWaveForm','suma_sg','t_total','samplingInterval','x')