% c_FixFsEx2   This example illustrates order tracking using fixed
%              sampling frequency on a signal with some tones passing
%              a resonance in an SDOF system. 

% This example is similar to the analysis behind figures 12.13 and 12.14 in Brandt

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
clc
close all

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

%--------------------------------------------------
% Define parameters
fs=2000;        		% Sampling frequency
T=30;					% Time duration
f0=10;					% Start freq.
f1=100;					% End freq.
t=(0:1/fs:T)';			% Time axis
Sq2=sqrt(2);			% Fundamental amplitude
% Define a tacho signal. Technically speaking, this signal does not contain
% 'pulses', as it is a sweeping sine. However, setting a trigger level it
% will be equivalent to a pulse train.
tacho=chirp(t,f0,t(end),f1);   
% Smoothing filter length (try changing this between 0 and 20 and see the
% result!)
FiltLength=10;          
% Compute an rpm profile
TrigLevel=0; Slope=1;
PPR=1;
[rpm,trpm]=tacho2rpm(tacho,fs,TrigLevel,Slope,PPR,fs,FiltLength);
% Plot the rpm vs. time
plot(t,rpm,LineType{1},'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('RPM','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

% Define a synthetic vibration signal with three instantaneouss frequencies: 
% The fundamental, and 2* and 3* this frequency
x=Sq2*tacho;                          % Tacho signal contains fundamental
x=x+0.5*Sq2*chirp(t,2*f0,t(end),2*f1);  % Add order 2
x=x+0.25*Sq2*chirp(t,3*f0,t(end),3*f1); % Add order 3
% Pass this signal through an SDOF system
fr=50; zr=0.02;
m=1; k=(2*pi*fr)^2*m;
c=2*zr*sqrt(k*m);
x=1e4*timefresp(x,fs,m,c,k,1,1,'d');

% Now produce a RPM map of this vibration signal, using the rpm-time
% profile from above
MinRpm=600; MaxRpm=6000;DeltaRpm=50;
Win=ahann(1024);
% Win=aflattop(1024);
[S,F,R] = rpmmap(x,fs,rpm,trpm,MinRpm,MaxRpm,DeltaRpm,Win);

% Next plot in color map format
figure
plotrpmmapc(S,F,R,4,1,'log')
ylim([0 500])

% Finally extract orders 1, 2, and 3 and plot
o=map2order(S,F,R,[1:3],'Hann5');
% o=map2order(S,F,R,[1:3],'peak');
figure
plot(R,o)
xlabel('RPM')
ylabel('Acceleration [m/s^2 RMS]')
axis([600 6000 0 3])
legend('Order 1','Order 2','Order 3')