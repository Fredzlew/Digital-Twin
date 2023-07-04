% f_SamplingPulse   This example illustrates the effect of bandwidth on
%                   representing a pulse in time domain

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

% This example produces a plot similar to Figure 3.4
% Set up parameters
N=1000;                  	% Number of samples to start with
fs=10000;                 	% Sampling freq. for simulation
tp=1e-3;                    % Half sine pulse time
% Use a toolbox command, makepulse, to create a half sine pulse
x=makepulse(N,fs,tp,'halfsine');
% Add zeros to the pulse
x=[zeros(N,1);x];
tx=makexaxis(x,1/fs);
% Now resample the high resolution pulse, which is equal to sampling the
% signal using an almost 'ideal' measurement system using a sampling
% frequency of 200 Hz (10000/50)
y=resample(x,1,50);
ty=makexaxis(y,1/(fs/50));
% Plot. Note the difference in time scales! Also note the drastical drop in
% amplitude, which is caused by the lowpass filtering removing most of the
% energy content of the original signal.
subplot(1,2,1)
plot(1000*tx,x,'k','LineWidth',LineWidth)
axis([99 102 -.05 .6])
title('Pulse sampled with 10 kHz sampling freq.','FontName',FontName,'FontSize',FontSize)
% legend('Original','Resampled')
xlabel('Time [ms]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(1,2,2)
plot(1000*ty,y,'k','LineWidth',LineWidth)
axis([50 150 -.005 .06])
xlabel('Time [ms]','FontName',FontName,'FontSize',FontSize)
title('Pulse sampled with 200 Hz sampling freq.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)

