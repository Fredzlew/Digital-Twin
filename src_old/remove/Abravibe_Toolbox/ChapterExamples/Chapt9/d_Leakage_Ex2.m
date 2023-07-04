% d_Leakage_Ex2     Example of leakage
%
% This example illustrates leakage in the DFT using rectangular (= no window), 
% Hanning, and flattop windows. 
% The example shows the broadening of the peak caused by the window; the
% smaller the amplitude uncertainty the window causes, the wider the peak.
% The plots are similar to those in Fig. 9.11

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
FontSize=11;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

%--------------------------------------------------
% Define a sine, periodic in time window, 32 periods
N=128; fs=1;
n=(0:1:N-1)';
y=sin(32*2*pi*n/N);
Y1=sqrt(2)*alinspec(y,fs,boxcar(N));      % Flattop
Y2=sqrt(2)*alinspec(y,fs,ahann(N));			% Hanning
Y3=sqrt(2)*alinspec(y,fs,aflattop(N));      % Flattop
k=(-N/4:N/4)';
figure
plot(n,y)
title('Time signal with 32 periods in time window')
figure
plot(k,Y1,'-ko',k,Y2,'-k+',k,Y3,'-k*')
xlabel('k','FontName',FontName,'FontSize',FontSize)
grid
axis([-8 8 0 1.1])
legend('Uniform','Hanning','Flattop')
title('DFT when sine is periodic in time window','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize-1)
% Now produce a similar sine, but with a half period more in the time
% window, which is identical to the fact its frequency is mid way between
% two frequency lines in the DFT
figure
y=sin(32.5*2*pi*n/N);
plot(n,y)
title('Time signal with 32.5 periods in time window')  % Bug fixed 26.9.2012! Said 32 before...
% The apparent "amplitude modulation" you see in the time plot, is simply
% caused by the few time samples, which means each period is sampled at a
% different phase angle of the sine. The information is intact, as is seen
% in the coming DFT plot:
figure
Y1=sqrt(2)*alinspec(y,fs,boxcar(N));
Y2=sqrt(2)*alinspec(y,fs,ahann(N));			% Hanning
Y3=sqrt(2)*alinspec(y,fs,aflattop(N));
% plot(k,Y1,'k',k,Y1,'ko',k,Y2,'k',k,Y2,'k+',k,Y3,'k',k,Y3,'k*')
plot(k,Y1,'-ko',k,Y2,'-k+',k,Y3,'-k*')
xlabel('k','FontName',FontName,'FontSize',FontSize)
grid
axis([-8 8 0 1.1])
legend('Uniform','Hanning','Flattop')
title('DFT when sine is NOT periodic in time window','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize-1)
