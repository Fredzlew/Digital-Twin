% b_DFT_Example2    Example of DFT
%
% This example illustrates the periodic repetition of the DFT. 
% This example is a CORRECTED version of Figure 9.4 in Brandt,
% where the figure was wrong.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% 2012-01-14 modified to work with Octave

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
% Original sequence
n=(0:15)';
y=cos(4*pi*n/16);
Y=fft(y)/16;
% mirroring to the left of n=0
n1=-flipud(n(1:end))-1;
% Original n
n2=n;
% n from 16 to 31
n3=n+16;
Y=real(Y);          % Imaginary part is zero anyway, so only plot real part
stem(n1(1:8),Y(1:8),'MarkerSize',4,'LineWidth',LineWidth)
axis([-16 32 0 10])
hold on
stem(n1(9:end),Y(9:end),'filled')
stem(n2(1:8),Y(1:8),'filled')
stem(n2(9:end),Y(9:end),'MarkerSize',4,'LineWidth',LineWidth)
stem(n3,Y,'MarkerSize',4,'LineWidth',LineWidth)
xlabel('Frequency number, k','FontName',FontName,'FontSize',FontSize)
ylabel('Real[X(k)]/N','FontName',FontName,'FontSize',FontSize)
axis([-16 32 -.1 1])
set(gca,'XTick',[-16 -8 0 8 16 24 32])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

fprintf('Please note that Figure 9.4, page 185 in Brandt is incorrect. Should look like the figure here!\n')
