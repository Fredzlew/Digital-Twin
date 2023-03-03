% e_EnvSpecEx   Example of calculation of envelope spectrum on an acceleration
%               signal measured on a milling machine during operation


% This example is similar to that used for Figure 16.10 in Brandt

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

% Load signal and compute envelope spectrum 
FileName='..\data\MillVib';
load(FileName);
[E,f]=aenvspec(x,fs,16*1024,1543,200);
figure;
plot(f,E,LineType{1},'LineWidth',LineWidth)
xlim([0 200])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Envelop Spectrum [-]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
