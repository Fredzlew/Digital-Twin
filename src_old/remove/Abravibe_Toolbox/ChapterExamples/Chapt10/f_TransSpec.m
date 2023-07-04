% f_TransSpec   Transient spectrum of a half sine pulse.
%

% This example is similar to the plot in Fig. 10.9 in Brandt.

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
% Half sine pulse transient spectrum
fs=1000;
N=1024;
tp=10e-3;           % Pulse duration
A=100;
p=makepulse(N,fs,tp,'gaussian');
% Scale pulse to max peak of A
p=A*p/max(p);
% Compute transient spectrum
[T,f]=atranspec(p,fs);
plot(f,T,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Transient spectrum [Ns]','FontName',FontName,'FontSize',FontSize)
title('Transient spectrum of half sine pulse','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
