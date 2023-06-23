% a_SRSEx     Calculate shock response spectrum (SRS) of a half sine pulse
% 

% This example is similar to that used for Figure 16.2 in Brandt

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
close all
clc

%--------------------------------------------------
FontSize=11;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------
% Define parameters
Nf=100;
Q=10;
d=11e-3;
fs=1e4;
N=fs;
fmin=1; fmax=2000;
% Calculate a half-sine pulse using the command makepulse
p=makepulse(N,fs,d,'halfsine');
p=100*p/max(p);         % Normalize to peak value of 100 g
% Compute SRS assuming the pulse is acceleration measured in g's
[SRS,f]=asrs(p,fs,fmin,fmax,Nf,Q);
% Plot result
figure;
loglog(f,SRS,LineType{1},'LineWidth',LineWidth);
xlim([fmin fmax])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Shock Response Spectrum, Q=10 [g]','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

