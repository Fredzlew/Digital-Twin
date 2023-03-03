% i_PSD_MidexSignal     Compute PSD and cumulated PSD of a random signal from a
%                       3DOF system mixed with a sine.

% This example is similar to the plot in Fig. 10.20 in Brandt.

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
% PSD of 3DOF system
FileName='..\data\3dofmkz.mat';
if exist(FileName,'file') ~= 2
    fprintf('You first must run example i_3DOF_MKz in folder ''chapter6''\n')
else
    load(FileName)
    t=(0:1/fs:(length(y)-1)/fs)';
    y=y+sqrt(2)*sin(2*pi*80*t);
    [Gyy,f]=apsdw(y,fs,2*N);
    subplot(2,1,1)
    semilogy(f,Gyy,LineType{1},'LineWidth',LineWidth)
    xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
    ylabel('Acceleration PSD [(m/s^2)^2/Hz]','FontName',FontName,'FontSize',FontSize)
    axis([0 400 1e-5 200])
    grid
    set(gca,'YTick',[1e-4 1e-2 1 1e2])
    set(gca,'FontName',FontName,'FontSize',FontSize)
    subplot(2,1,2)
    C=cumsum(Gyy)*f(2);
    semilogy(f,C,LineType{1},'LineWidth',LineWidth)
    xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
    ylabel('Cumulated PSD [(m/s^2)^2]','FontName',FontName,'FontSize',FontSize)
    axis([0 400 1e-5 1000])
    grid
    set(gca,'YTick',[1e-4 1e-2 1 1e2])
    set(gca,'FontName',FontName,'FontSize',FontSize)
end