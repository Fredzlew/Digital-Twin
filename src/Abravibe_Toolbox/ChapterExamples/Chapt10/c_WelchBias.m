% c_WelchBias      example of bias in a PSD estimated using Welch's method
%                  on a displacement signal from an SDOF system excited by
%                  bandlimited white noise.
%

% This example is similar to the plot in Fig. 10.4 in Brandt.

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
% Set up a mechanical system
fr=10;             % Integer number!
zr=0.05;
m=1;
k=(2*pi*fr)^2*m;    % Resonance will be fr Hz
c=zr*2*sqrt(m*k);   % Relative damping will be zr
DType='d';          % 'd', 'v', 'a' for displ. velocity, acc.
%                   % Note: This example works best if you keep
%                   displacement as the unit of output data. Also, the axes
%                   texts do not adapt if you change this parameter.

% Define parameters for FFT analysis
fs=128;
N=256;
M=500;              % Number independent averages for Welch
% Compute a random input force to the mechanical system (rms =1 N)
F=1000*randn(M*N,1);
% ...and compute the DType output of the SDOF system
x=timefresp(F,fs,m,c,k,1,1,DType);

% Compute true PSD which is Pxx=PFF*abs(H)^2
PFF=1e6/(fs/2);       % Constant PSD, RMS=1 kN
f=(0:fs/N:fs/2)';   % Length N/2+1
H=mck2frf(f,m,c,k,1,1,DType);
Pxx=PFF*abs(H).^2;
% Compute Welch PSD
[Pxx_w,f]=apsdw(x,fs,N);
% Plot data
semilogy(f,Pxx_w,LineType{1},f,Pxx,LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Displacement PSD [m^2/Hz]','FontName',FontName,'FontSize',FontSize)
grid
xlim([fr-5 fr+5])
legend('Welch estimate','True PSD')
set(gca,'FontName',FontName,'FontSize',FontSize)
