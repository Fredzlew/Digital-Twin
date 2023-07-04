% a_DFT_Ex1     DFT example
%
% This example shows the DFT result of a cosine with two periods within the
% time window, in real and imaginary parts, respectively. Since a cosine is
% an even function, only the real part of the DFT contains nonzero data.
% This example is similar to the results in Figure 9.3 in Brandt.

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
% Create a cosine with 2 periods, N=16
n=(0:15)';
y=cos(4*pi*n/16);
% Spectrum
Y=fft(y)/16;
% Plot data
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
subplot(3,1,1)
stem(n,y,'fill','MarkerSize',4,'LineWidth',LineWidth)
hold on
n=[n;16];
plot(n,cos(4*pi*n/16),LineType{1},'LineWidth',LineWidth)
xlabel('Sample number, n','FontName',FontName,'FontSize',FontSize)
ylabel('x(n)','FontName',FontName,'FontSize',FontSize)
axis([0 16 -2 2])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,1,2)
n=n(1:end-1);
stem(n,real(Y),'fill','MarkerSize',4,'LineWidth',LineWidth)
xlabel('Frequency number, k','FontName',FontName,'FontSize',FontSize)
ylabel('Real[X(k)]/N','FontName',FontName,'FontSize',FontSize)
axis([0 16 -1 1])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,1,3)
stem(n,imag(Y),'fill','MarkerSize',4,'LineWidth',LineWidth)
axis([0 16 -1 1])
xlabel('Frequency number, k','FontName',FontName,'FontSize',FontSize)
ylabel('Imag[X(k)]/N','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
