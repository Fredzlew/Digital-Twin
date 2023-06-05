% i_FanSpectra  Compare two linear spectra with too coarse and sufficient
%               frequency resolution. It illustrates the need to increase
%               the blocksize (decrease frequency increment) until a clear
%               line spectrum is visible, in the case of (relatively clean)
%               periodic signals.
%

% This example is similar to the plot in Fig. 10.16 in Brandt.

% The data in this example come from measurements on a small fan with some
% imbalance causing vibrations with particularly first order.

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
% Spectra of fan
FileName1='..\data\fan1.mat';
FileName2='..\data\fan3.mat';
load(FileName1);
X1=Data;
f1=makexaxis(X1,Header.xIncrement);
load(FileName2);
X2=Data;
f2=makexaxis(X2,Header.xIncrement);
% Plot data
subplot(1,2,1)
plot(f1,X1,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Linear Spectrum Acc. [m/s^2 RMS]','FontName',FontName,'FontSize',FontSize)
grid
xlim([0 500])
title('Too coarse frequency increment','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(1,2,2)
plot(f2,X2,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Linear Spectrum Acc. [m/s^2 RMS]','FontName',FontName,'FontSize',FontSize)
grid
xlim([0 500])
title('Sufficient frequency increment','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
