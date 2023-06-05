% k_SmPerPSD    Compute PSD on data from a simulation of a 3DOF system,
%               using the smoothed periodogram method, wiht logarithmic
%               frequency resolution.

% This example is similar to the plot in Fig. 10.19 in Brandt.

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
    % Use the smoothed periodogram for PSD estimation
    fmin=40;
    fmax=round(fs/2.5);
    NFreqs=400;
    Lsmin=100;
    [Gyysp,fsp,Nw]=apsdsp(y,fs,fmin,fmax,NFreqs,Lsmin);
    % Plot results overlaid by true PSD
    semilogy(fsp,Gyysp,LineType{1},f,Gyyt,LineType{4},'LineWidth',LineWidth)
    xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
    ylabel('Acceleration PSD [(m/s^2)^2/Hz]','FontName',FontName,'FontSize',FontSize)
    axis([40 400 1e-5 100])
    grid
    set(gca,'YTick',[1e-4 1e-2 1 1e2])
    set(gca,'FontName',FontName,'FontSize',FontSize)
end