% a_FRFandIR        Calculate frequency response and impulse response of 
%                   SDOF system

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

% Define mechanical system
m=1;
c=100;
k=1e6;
f0=1/2/pi*sqrt(k/m);        % f0=159.2 Hz
z=c/2/sqrt(m*k);            % z=5 %

% Create a suitable frequency axis for computation of the frequency
% response
f=(0:500/2049:500-1/2049)';
% Now calculate the FRF and its IR
Hd=mck2frf(f,m,c,k,1,1,'d');
[h5,t]=frf2ir(Hd,f);

% Plot frequency response
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
semilogy(f,abs(Hd),'k','LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Frequency Response [m/N]','FontName',FontName,'FontSize',FontSize)
% axis([0 .5 -1e-3 1e-3])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

% Plot impulse response
% Cut impulse responses at 1 sec
idx=min(find(t >= 1));
t=t(1:idx);
h5=h5(1:idx);
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
plot(t,h5,'k','LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Impulse Response [m/Ns]','FontName',FontName,'FontSize',FontSize)
% axis([0 .5 -1e-3 1e-3])
grid
axis tight
set(gca,'FontName',FontName,'FontSize',FontSize)