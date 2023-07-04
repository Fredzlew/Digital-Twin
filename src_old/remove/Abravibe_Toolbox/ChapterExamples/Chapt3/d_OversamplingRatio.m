% d_OversamplingRatio       This example shows the effect of the oversampling 
%                           ratio on the appearance of a signal


% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.


close all
clear

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------

%================================================================
% Create a random signal with 50 samples. Oversampling ratio is exactly 2.
x=randn(50,1);
% Resample to an oversampling ratio of 2.5
x=resample(x,4,5);
x=resample(x,5,4);      % This now has oversampling ratio 2.5, comparable to 
                        % most measured signals from noise and vibration
                        % instruments
% Resample x up to 25 times oversampling
fs=1000;                % To make the example with a reality feeling, we 
                        % assume x has a sampling frequency in Hz
x2=resample(x,10,1);
t=makexaxis(x,1/fs);
t2=makexaxis(x2,1/(10*fs));
% Plot the upsampled signal to show it is smoothly varying
subplot(1,2,1)
plot(t,x,'k','LineWidth',LineWidth)
axis([0 0.05 -2.5 2.5])
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
title('Signal with 2.5 times oversampling','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(1,2,2)
% ... and compare to the signal with less oversampling. Both signals
% fulfill the sampling theorem, but the time domain appearance of the low
% oversampling ratio of 2.5 does not give a good description of what the
% signal actually looks like.
plot(t2,x2,'k','LineWidth',LineWidth)
axis([0 0.05 -2.5 2.5])
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
title('Signal with 25 times oversampling','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
