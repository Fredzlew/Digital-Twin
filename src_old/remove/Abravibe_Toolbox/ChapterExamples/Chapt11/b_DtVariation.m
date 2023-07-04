% b_DtVariation     This example illustrates the effect of limited dynamic
%                   range caused by scattering in the sampling moments, that
%                   is, the time between samples is not constant.

% This example is similar to the plot in Fig. 11.4 in Brandt.

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
% Set the standard deviation of the error in delta_t, relative to delta_t
% itself
target_std_for_error=1e-4;
% Set simulation parameters
fs=1024;
f0=300;
N=32*1024;
t=(0:1/fs:N/fs-1/fs)';	% 100 Hz sine, sampled with 1 kHz, N=1024
% Define an error vector with errors in seconds
err=target_std_for_error/fs*(randn(size(t)));
y=sqrt(2)*sin(2*f0*pi*t);	% "True" sampling at equidistant times
% Now make another time axis with error in the sampling instants
t_err=t.*1+err;
% ...and compute the sine at these erroneous time instants
y_err=sqrt(2)*sin(2*f0*pi*t_err);
% Compute a spectrum of the results, both true and erroneous
[A,f]=alinspec(y,fs,boxcar(N),1,0);						% Linear spectrum on perfect signal
[A_err,f]=alinspec(y_err,fs,boxcar(N),1,0);             % Linear spectrum on error signal
% Plot the results in dB scale: The dynamic range reduces from the approx.
% 300 dB caused by MATLABs double precision, to approx. 100 dB
figure
plot(f,[db20(A) db20(A_err)])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Lin. Spectrum [dB rel. 1 V rms]','FontName',FontName,'FontSize',FontSize)
title('Reduction in dynamic range due to erroneous sampling instants','FontName',FontName,'FontSize',FontSize)
grid
axis([0 512 -350 10])
set(gca,'FontName',FontName,'FontSize',FontSize)
