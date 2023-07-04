% c2_Resampling       Illustrates the interpolation formula

% This example produces a plot similar to Figure 3.6 in the
% book. It illustrates the process of resampling a signal to a higher
% sampling frequency.

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

%--------------------------------------------------
% Set up
% Note that because we use random noise (randn command) each run of this
% script will produce a unique realization. Sometimes the result may look
% better than other times.
N=1000;             % Number of samples to start with
fs=1000;           	% Sampling freq. for simulation
x=randn(N,1);      	% Gaussian noise, oversampling ratio is 2
x=resample(x,5,1);  % Resample to 10 times oversampling
x25=x(1:4:end);     % Oversampling ratio of x25 is 2.5
xr=resample(x25,4,1);	% oversampling ratio of xr is 10
% Plot 'original' data x, and the resampled xr
plot(1:50,x(1:50),'ok',1:50,xr(1:50),'+k','LineWidth',LineWidth)
legend('Original','Resampled',3)
xlabel('Sample Number','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
