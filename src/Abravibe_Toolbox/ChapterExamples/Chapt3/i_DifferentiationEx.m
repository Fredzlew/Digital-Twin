% i_DifferentiationEx.m Example of differentiation of time velocity into
% acceleration
%
% For reference, see the book section 3.4.3.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% NOTE for Octave Users! There currently seems to be a problem with the
% remez command in Octave, which means that timediff does not work as
% expected.

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

% Create a sine to be differentiated
t=makexaxis([1:N]',1/fs);
x=sin(2*pi*f0*t);

% Compute the true integral of x, which is a sine
yt=(2*pi*f0)*cos(2*pi*f0*t);
tt=makexaxis(yt,1/fs);

% Now use the ABRAVIBE command TIMEDIFF to differentiate it. See inside 
% TIMEDIFF to understand how it does the differentiation. 
Type='remez';               
y = timediff(x,fs,Type);
% Note! timediff changes the length of y
t=t(1:length(y));
figure
plot(tt,yt,t,y);
legend('True sine','Integral of cosine')
title(['Result of timediff with option '' ' upper(Type) ' '' '])

Type='maxflat';               
y = timediff(x,fs,Type);
% Note! timediff changes the length of y
t=makexaxis([1:N]',1/fs);   % Original time axis again (was truncated for type remez above)
t=t(1:length(y));
figure
plot(tt,yt,t,y);
legend('True sine','Integral of cosine')
title(['Result of timediff with option '' ' upper(Type) ' '' '])

