% e_SynchFsEx2  This example illustrates order tracking using
%               synchronuous sampling frequency, including phase. 
%               An rpm-time profile is created, then the time signal is 
%               resampled, after which a color map plot, and finally 
%               orders vs. RPM are extracted in ampitude/phase format.

% This example is similar to d_SynchEx2 but also shows phase of orders

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
fs=2*1024;        		% Sampling frequency
T=90;					% Time duration
f0=5;					% Start freq.
f1=100;					% End freq.
t=(0:1/fs:T)';			% Time axis
Sq2=sqrt(2);			% Fundamental amplitude
% Define a tacho signal. Technically speaking, this signal does not contain
% 'pulses', as it is a sweeping sine. However, setting a trigger level it
% will be equivalent to a pulse train.
tacho=chirp(t,f0,t(end),f1);   
CutLevel=0.6;           % Add some harmonics by limiting amplitude
tacho(tacho>CutLevel)=CutLevel;
tacho(tacho<-CutLevel)=-CutLevel;
% Smoothing filter length (try changing this between 0 and 20 and see the
% result!)
FiltLength=20;          
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

% Define a synthetic vibration signal with three instantaneous frequencies: 
% The fundamental, and 2* and 3* this frequency
x=Sq2*tacho;                            % Tacho signal contains fundamental
x=x+x.^3;
x(x>CutLevel*max(x))=CutLevel*max(x);
x(x<-CutLevel*max(x))=-CutLevel*max(x);% Pass this signal through an SDOF system
fr=60; zr=0.02;
m=1; k=(2*pi*fr)^2*m;
c=2*zr*sqrt(k*m);
x=1e4*timefresp(x,fs,m,c,k,1,1,'d');

% This vibration signal is now resampled using the rpm-time
% profile from above
MaxOrd=32;
[xs,ts] = synchsampr(x,fs,rpm,MaxOrd);
% Also resample tacho signal for phase ref:
[tachos,ts] = synchsampr(tacho,fs,rpm,MaxOrd);

% The resampled signal is now processed into an RPM map
MinRpm=600; MaxRpm=6000;DeltaRpm=10;
OrdRes=16;
WinStr='aflattop';
[S,F,R] = rpmmaps(xs,rpm,OrdRes,MaxOrd,MinRpm,MaxRpm,DeltaRpm,WinStr,tachos);

% Next plot in color map format
figure
plotrpmmapsc(S,F,R,3,1,'log')

% Finally extract orders 1, 2, and 3 and plot
% o=map2order(S,F,R,[1:3],'Hann5');
o=map2order(S,F,R,[1 2],'synch');
figure
subplot(2,1,1)
plot(R,abs(o))
ylabel('Acceleration [m/s^2 RMS]')
axis([600 6000 0 3])
legend('Order 1','Order 3')
subplot(2,1,2)
plot(R,angledeg(o+pi))
xlabel('RPM')
ylabel('Phase [Deg.]')