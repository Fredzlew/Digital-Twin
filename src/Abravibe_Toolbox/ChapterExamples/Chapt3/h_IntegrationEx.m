% h_IntegrationEx.m Example of integration of time acceleration into velocity
%
% For reference, see the book section 3.4.2.

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
% Define parameters
fs=1024;
N=20*fs;                  % Number of samples = 20 seconds
f0=17.1;                 % frequency of sine

% Create a cosine to be integrated
t=makexaxis([1:N]',1/fs);
x=cos(2*pi*f0*t);

% Compute the true integral of x, which is a sine
yt=1/(2*pi*f0)*sin(2*pi*f0*t);
tt=makexaxis(yt,1/fs);

% Now use the ABRAVIBE command TIMEINT to integrate it. See inside TIMEINT
% to understand how it does the integration. 
Type='simple';                 % Change to 'hpfilter', or 'simple' to 
                            % investigate differences 
y = timeint(x,fs,Type);
figure
plot(tt,yt,t,y);
legend('True sine','Integral of cosine')
title(['Result of timeint with option '' ' upper(Type) ' '' '])
y=y(1:end-N);             % Remove last block that may contain transients
L=length(y);

Type='hpfilter';            
Par=1;
y = timeint(x,fs,Type,Par);
figure
plot(tt,yt,t,y);
legend('True sine','Integral of cosine')
title(['Result of timeint with option '' ' upper(Type) ' '' '])
y=y(1:end-N);             % Remove last block that may contain transients
L=length(y);

Type='iir';            
Par=1;
y = timeint(x,fs,Type,Par);
% Note! option 'iir' removes some samples from the output signal due to
% time delays, to make the output y time synchronous with x
L=length(y);
t=t(1:L);
figure
plot(tt,yt,t,y);
legend('True sine','Integral of cosine')
title(['Result of timeint with option '' ' upper(Type) ' '' '])
y=y(1:end-N);             % Remove last block that may contain transients
L=length(y);
