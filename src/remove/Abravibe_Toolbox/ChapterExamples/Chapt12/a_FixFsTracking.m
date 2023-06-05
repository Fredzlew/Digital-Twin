% a_FixFsTracking   This example illustrates order tracking using fixed
%                   sampling frequency. An rpm-time profile is created, then
%                   a waterfall diagram of the RPM map, followed by a color
%                   map plot, and finally orders vs. RPM are extracted.

% This example is similar to example 12.3.1 in Brandt

% After running this example you may want to try the following changes:
% 1. Try using a flattop window with rpmmap, followed by option 'peak' in
%    the map2order command
% 2. Using the default parameters again, try using different block sizes,
%    512, 2048, and 4096 in the rpmmap command (by changing the length of
%    Win). Try to understand the results. 

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

% Define a synthetic vibration signal with three instantaneous frequencies: 
% The fundamental, and 2* and 3* this frequency
x=Sq2*tacho;                          % Tacho signal contains fundamental
x=x+0.5*Sq2*chirp(t,2*f0,t(end),2*f1);  % Add order 2
x=x+0.25*Sq2*chirp(t,3*f0,t(end),3*f1); % Add order 3
% Now produce a RPM map of this vibration signal, using the rpm-time
% profile from above
MinRpm=600; MaxRpm=6000;DeltaRpm=50;
Win=ahann(1024);
% Win=aflattop(1024);
[S,F,R] = rpmmap(x,fs,rpm,trpm,MinRpm,MaxRpm,DeltaRpm,Win);

% First plot in waterfall format
figure
plotrpmmapwf(S,F,R)
axis tight
xlim([0 500])
drawnow

% It may be useful to plot a single spectrum at a particular rpm:
idx=min(find(R >= 3000));    % Corrected 2013-09-25 make sure contains value
figure
plot(F,S(:,idx));
xlabel('Frequency [Hz]')
ylabel('Acceleration [m/s^2 RMS]')
title(['Spectrum of RPM map for ' num2str(R(idx)) ' RPM'])

% Next plot in color map format
figure
plotrpmmapc(S,F,R)

% Finally extract orders 1, 2, and 3 and plot
o=map2order(S,F,R,[1:3],'Hann5');
% o=map2order(S,F,R,[1:3],'peak');
figure
plot(R,o)
xlabel('RPM')
ylabel('Acceleration [m/s^2 RMS]')
axis([600 6000 0 1.1])