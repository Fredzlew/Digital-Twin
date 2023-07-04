% b_Beating    Example of beating in the response of an SDOF system
%
% Calculate frequency response and impulse response of SDOF system excited
% by a sine force with a frequency near the natural frequency of the SDOF
% system. The result is a beating phenomenon.

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

% Define mechanical system with low damping
m=1;
c=4;
k=1e6;
f0=1/2/pi*sqrt(k/m);        % f0=159.2 Hz
z=c/2/sqrt(m*k);            % z=0.2 %

%=========================================================================
% SDOF forced response with  beating. To simulate this we use the forced
% response time domain function TIMEFRESP
f01=155;            % Sine frequency close to f0
fs=2000;
t=(0:1/fs:1)';
% Create input signals with frequencies f01 and f02
x1=10*sin(2*pi*f01*t);
% Sum the forced response due to each input
y=timefresp(x1,fs,m,c,k,1,1,'d');
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gcf, 'PaperSize', [12 6])
plot(t,y,LineType{1},'LineWidth',LineWidth);
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Displacement [m]','FontName',FontName,'FontSize',FontSize)
title('Forced response of SDOF, f0=159.2 Hz, \zeta=0.2%, force freq.=155 Hz','FontName',FontName,'FontSize',FontSize)
grid

