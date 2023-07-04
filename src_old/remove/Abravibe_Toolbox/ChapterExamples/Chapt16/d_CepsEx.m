% d_CepsEx      Example of calculation of power cepstrum on an acceleration
%               signal measured on a milling machine during operation


% This example is similar to that used for Figure 16.8 in Brandt

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

% Load signal and compute linear spectrum and power cepstrum
FileName='..\data\MillVib';
load(FileName);
t=makexaxis(x,1/fs);
[X,f]=alinspec(x,fs,ahann(32*1024),1);
[cp,tc]=apceps(X.^2,f);
% Plot results
figure;
subplot(2,2,1)
plot(t,x,LineType{1},'LineWidth',LineWidth)
axis([0 0.05 -6 6])
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Acceleration [m/s^2]','FontName',FontName,'FontSize',FontSize)
title('Time signal','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
plot(f,X,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Lin. Spec [m/s^2 RMS]','FontName',FontName,'FontSize',FontSize)
axis([0 3000 0 0.5])
title('Spectrum','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
semilogy(f,X,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Lin. Spec [m/s^2 RMS]','FontName',FontName,'FontSize',FontSize)
axis([0 3000 1e-6 1])
title('Spectrum','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
plot(tc,abs(cp),LineType{1},'LineWidth',LineWidth)
xlim([0 .2])
xlabel('Quefrency [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Cepstrum [-]','FontName',FontName,'FontSize',FontSize)
title('Power Cepstrum','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
