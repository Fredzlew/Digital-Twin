% h_LinSpecVsPSD2 Compare a linear spectrum with a PSD of a random signal
%

% This example is similar to the plot in Fig. 10.14 in Brandt.

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
% Linspec vs PSD of sine
fs=6400;
N1=2*512;
N2=512;
f0=200;
t=(0:1/fs:(N1-1)/fs)';
y=randn(100*N1,1);
[X1,f1]=alinspec(y,fs,ahann(N1),100);
[X2,f2]=alinspec(y,fs,ahann(N2),100);
[P1,f1]=apsdw(y,fs,N1,0);
[P2,f2]=apsdw(y,fs,N2,0);
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 6]);
set(gcf, 'PaperSize', [12 6])
subplot(1,2,1)
plot(f1,X1,LineType{1},f2,X2,LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Linear Spectrum [V RMS]','FontName',FontName,'FontSize',FontSize)
grid
axis([0 1000 0 .12])
title('a)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(1,2,2)
plot(f1,P1,LineType{1},f2,P2,LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Spectral density [V^2/Hz]','FontName',FontName,'FontSize',FontSize)
grid
%xlim([150 250])
axis([0 1000 0 5e-4])
title('b)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)

% Now, for interpretation:
% In the linear spectra, you should compare the peaks, as they are
% interpreted as a sine with that RMS level:
fprintf('From the time signal we calculate the RMS level to %6.2f V\n',std(y))
R1=sqrt(sum(X1.^2)/winenbw(ahann(64)));
R2=sqrt(sum(X2.^2)/winenbw(ahann(64)));
fprintf('The RMS levels from the two linear spectra are: %6.2f, and %6.2f, V respectively\n',R1,R2);
% From the PSDs, we need to calculate the RMS level as the area under the
% PSDs between 0 and fs/2
df1=f1(2)-f1(1);
df2=f2(2)-f2(1);
R3=sqrt(df1*sum(P1));
R4=sqrt(df2*sum(P2));
fprintf('The RMS levels from the two PSDs are: %6.2f, and %6.2f, V respectively\n',R3,R4);
