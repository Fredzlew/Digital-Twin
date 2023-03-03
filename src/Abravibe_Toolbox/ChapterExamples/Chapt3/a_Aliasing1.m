% Aliasing 1. This example shows a simple illustration of aliasing
% The example samples two cosines with frequencies 90 and 110 Hz,
% respectively, with a sampling frequency of 200 Hz. The example
% illustrates the concept of aliasing, which says that, in this case with a
% Nyquist frequency of 100 Hz, a cosine of 110 Hz will alias as a cosine
% with 90 Hz. 


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
% This example illustrates the effect of 'sampling' a cosine signal with a
% sampling frequency (or time increment) which does not satisfy the
% sampling theorem.

% First we assume a frequency of the cosine. Let us say it is 90 Hz and we
% describe it well, with say 100 times oversampling to have as the 'true'
% signal
f1=90;
fs=100*f1;
t=(0:1/fs:10/f1);
x=cos(2*pi*f1*t);

% Next, we assume we sample this signal with a sampling frequency which is
% slightly more than needed. Say fs1=200 Hz which means the Nyquist
% frequency is 100 Hz
fs1=200;
t1=(0:1/fs1:10/f1);
y1=cos(2*pi*f1*t1);     % This is the 90 Hz cosine sampled with 200 Hz

% The aliasing phenomenon is now illustrated by instead sampling a 110 Hz
% cosine with the same 200 Hz sampling frequency. Note that the 90 Hz
% signal is as much below the Nyquist frequency, as the 110 Hz is above.
% This is what happens:
f2=110;
y2=cos(2*pi*f2*t1);     % The 110 Hz cosine sampled with 200 Hz

% Now we plot the three signals. The 'true', the correctly sampled 90 Hz,
% and the incorrectly sampled 110 Hz signals.
% We use green plus signs (+) for the correct signal, and red rings (o) for
% the incorrectly sampled signal
plot(t,x,'k',t1,y1,'g+',t1,y2,'ro','LineWidth',LineWidth)
legend('True','90 Hz','110 Hz','FontName',FontName,'FontSize',FontSize)
title('Note that the 90 and 110 Hz cosines have identical values at the samples','FontName',FontName,'FontSize',FontSize)