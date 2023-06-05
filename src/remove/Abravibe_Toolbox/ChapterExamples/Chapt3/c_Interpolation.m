% c_Interpolation       Illustrates the interpolation formula

% This example produces a plot similar to Figure 3.5 in the
% book. It illustrates the process of interpolation of a signal as given
% by Eq. (3.2).
% The created sample is, in this example, a relatively crude approximation,
% as only 65 samples are used to create it (some outside the plot created
% by this script! remove the axis(...) command to see the entire data used).

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
% Set up example
% This example uses a nominal sampling frequency of one Hz for
% simplicity. It shows how a new sample at time t0 can be computed by the
% interpolation formula in Eq. (3.2), from a number of samples of a signal.
t0=0.25;                    % Time of new sample
% Produce a discrete axis with 65 samples, numbered -32 to +32,
% corresponding to the sampling instances we simulate
n=[-32:32]';
% Produce a fine time axis to plot 'continuous' functions
t=(-32:.01:32)';
% Produce the true signal, a sine with 0.25 Hz frequency
x=0.5*sin(2*pi*0.25*t);
% Sample this true signal with fs=1 Hz
xn=0.5*sin(2*pi*0.25*n);
% Now calculate the time shifted sinc function in 'continuous' time, i.e.
% on time axis t:
s=sinc(t-t0);
% and next 'sample' this at the sampling instances in variable n
sn=sinc(n-t0);
%==========================================================================
newx=sum(sn.*xn);           % This calculates the new interpolated sample
%                           % as by Eq. (3.2)
%==========================================================================
% Plot 'original' data x, and the resampled xr
p=plot(n,xn,'o',n,sn,'^',t0,newx,'*',t,x,'-',t,s,'-');
set(p,'markerfacecolor','k','MarkerSize',5)
axis([-8 8 -1.1 1.1])
grid
legend('Signal Sample','Calculated Sinc Sample','New Sample','Location','SouthWest')
xlabel('Sample Number','FontName',FontName,'FontSize',FontSize)
title('Illustration of creating a new sample at x=0.25, from a number of other samples','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
