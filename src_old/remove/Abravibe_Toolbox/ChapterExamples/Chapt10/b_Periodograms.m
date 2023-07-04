% b_Periodograms        Example of periodograms of random data, and random data
%                       with sine tones hidden
%
% This example illustrates periodograms of two different lengths, in two cases:
% 1) The data are random noise, and 
% 2) The data are random noise plus periodic components
% You should note how the periodic components stand out in the long
% blocksize periodogram (lower right plot) compared to the shorter
% blocksize. You should also note how the pure random noise data does not
% improve with the longer blocksize, rather the opposite: The periodogram
% is thus not a good estimate of a PSD but very useful for periodic signals 
% hidden in noise.
%
% This example is similar to the plot in Fig. 10.2 in Brandt.

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
% Set parameters of the two blocksizes
fs=1000;
N1=512;
N2=8*N1;
% Create time axes for periodic signals
t1=(0:1/fs:(N1-1)/fs)';
t2=(0:1/fs:(N2-1)/fs)';
% Set frequency of square wave
f1=fs/30;
y1=2*square(2*pi*f1*t1);
y2=2*square(2*pi*f1*t2);
% Create two random sequences
x1=randn(N1,1);
x2=randn(N2,1);
% Create periodograms of both signals
P1=abs(fft(x1)).^2/N1;
P1=P1(1:N1/2+1);
P2=abs(fft(x2)).^2/N2;
P2=P2(1:N2/2+1);
f1=(0:fs/N1:fs/2);
f2=(0:fs/N2:fs/2);
% Add the periodic components to each random signal and compute
% periodograms
x3=x1+y1;
P3=abs(fft(x3)).^2/N1;
P3=P3(1:N1/2+1);
x4=x2+y2;
P4=abs(fft(x4)).^2/N2;
P4=P4(1:N2/2+1);
% Plot data
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
subplot(2,2,1)
semilogy(f1,P1,'k')
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel(['Periodogram length N= ',int2str(N1)],'FontName',FontName,'FontSize',FontSize)
grid
title('Random data only','FontName',FontName,'FontSize',FontSize)
axis([0 500 1e-5 10])
set(gca,'YTick',[1e-5 1e-3 1e-1 1e1])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
semilogy(f2,P2,'k')
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel(['Periodogram length N= ',int2str(N2)],'FontName',FontName,'FontSize',FontSize)
grid
axis([0 500 1e-5 10])
set(gca,'YTick',[1e-5 1e-3 1e-1 1e1])
title('Random data only','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
semilogy(f1,P3,'k')
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel(['Periodogram length N= ',int2str(N1)],'FontName',FontName,'FontSize',FontSize)
grid
axis([0 500 1e-2 1e5])
set(gca,'YTick',[1e-1 1e1 1e3 1e5])
title('Random data plus periodic components','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
semilogy(f2,P4,'k')
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel(['Periodogram length N= ',int2str(N2)],'FontName',FontName,'FontSize',FontSize)
grid
axis([0 500 1e-2 1e5])
set(gca,'YTick',[1e-1 1e1 1e3 1e5])
title('Random data plus periodic components','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
