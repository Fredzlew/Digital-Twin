% a_CumPSD  	Calculate cumulated mean square value of a PSD

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
% Start by loading a PSD from a measurement of acceleration on a plate
FileName='..\Data\PlatePSD';
load(FileName)
%================================================================
% PSD of plate vibration
subplot(2,1,1)
semilogy(f,Data,'k','LineWidth',LineWidth)
axis([0 2000 1e-3 2])
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Acceleration PSD [(m/s^2)^2/Hz]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)

%================================================================
% Figure 8.6: Cumulated spectrum
df=f(2)-f(1);
subplot(2,1,2)
C=cumsum(Data)*df;
semilogy(f,C,'k','LineWidth',LineWidth)
axis([0 2000 1e-2 200])
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Cumulated PSD [(m/s^2)^2]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
