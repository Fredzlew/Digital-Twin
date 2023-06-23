% c_Leakage_Ex1         Leakage example
%
% This example illustrates leakage in the DFT. 

% This example is similar to the results in Figure 9.5 in Brandt.

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
% Essentially same as Figure 9.5:
t=(0:1/512:.5-1/512)';
y=sin(2*pi*51*t);
Y=2*abs(fft(y))/256;
Y=Y(1:129);
f=(0:2:256)';
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
subplot(2,1,1)
plot(t,y,LineType{1},'LineWidth',LineWidth);
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('y=sin(2\pi51t)','FontName',FontName,'FontSize',FontSize)
grid
title('Leakage')
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,1,2)
plot(f,Y,LineType{1},'LineWidth',LineWidth);
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('2 \cdot |fft(y)/N|','FontName',FontName,'FontSize',FontSize)
grid
axis([0 256 -.05 1])
set(gca,'FontName',FontName,'FontSize',FontSize)
