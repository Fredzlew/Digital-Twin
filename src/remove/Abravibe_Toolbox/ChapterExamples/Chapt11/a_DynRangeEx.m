% a_DynRangeEx  This example illustrates the effect of limited dynamic
%               range caused by using too few bits in the ADC, that is,
%               setting too large full scale range.


% This example is similar to the plot in Fig. 11.3 in Brandt.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% 2012-01-15 Modified to work with Octave

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
fs=1000;            % Sampling frequency
N=2048;             % FFT block size
fsine=67;
t=(0:1/fs:(N-1)/fs)';   % Time axis
y=sqrt(2)*cos(2*pi*fsine*t);       % fsine Hz sine, rms value = 1
% Now let us simulate sampling of the sine in y, with 16 bits resolution
% and a full scale voltage set to 2 volts (little over the max amplitude of
% sqrt(2):
warning off
ScaleFactor=double(intmax('int16')/2);     % This makes 2V = maxint
ys1=double(int16(ScaleFactor*y));           % Truncate product to 16 bits
[Y,f]=alinspec(ys1/ScaleFactor,fs,ahann(N),1,0);    % Plot vs. original y scale
% Plot spectrum
hf=figure;
% subplot(1,2,1,'Parent',hf,'YScale','log','YMinorTick','off',...
%     'YMinorGrid','off',...
%     'FontSize',9,...
%     'FontName','Times New Roman')
subplot(1,2,1)
semilogy(f,Y,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Lin. Spectrum [V rms]','FontName',FontName,'FontSize',FontSize)
title('Full scale range 2 V','FontName',FontName,'FontSize',FontSize)
grid
axis([0 400 1e-7 5])
% set(gca,'YScale','log')
set(gca,'YTick',[1e-7 1e-5 1e-3 1e-1 1e0])
set(gca,'FontName',FontName,'FontSize',FontSize)
% Next, simulate measuring the same signal at full scale voltage of 20 V
ScaleFactor=double(intmax('int16')/20);     % This makes 20V = maxint
ys2=double(int16(ScaleFactor*y));
[Y,f]=alinspec(ys2,fs,ahann(N),1,0);
warning on
% Plot spectrum
% subplot(1,2,2,'Parent',hf,'YScale','log','YMinorTick','off',...
%     'YMinorGrid','off',...
%     'FontSize',9,...
%     'FontName','Times New Roman')
subplot(1,2,2)
semilogy(f,Y/ScaleFactor,LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Lin. Spectrum [V rms]','FontName',FontName,'FontSize',FontSize)
title('Full scale range 20 V','FontName',FontName,'FontSize',FontSize)
grid
axis([0 400 1e-7 5])
set(gca,'YTick',[1e-7 1e-5 1e-3 1e-1 1e0])
set(gca,'FontName',FontName,'FontSize',FontSize)
