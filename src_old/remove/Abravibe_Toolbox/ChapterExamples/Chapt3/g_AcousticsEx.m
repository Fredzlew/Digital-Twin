% g_AcousticsEx.m    Third-octave filtering and SLM example
%
% This example shows how to utilize some ABRAVIBE toolbox commands to
% understand and use 1/n octave filters and emulate a sound level meter 
% (SLM). For reference, see the book sections 3.3.4 - 3.3.6.

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
% Settings
fs=44100;               % Sampling frequency. Common audio frequency
n=3;                    % 1/n octave filter
fc=500;                 % Center frequency of 1/n octave filter
%--------------------------------------------------
% First example here will create coefficients for a third-octave filter for
% the 500 Hz third octave.
% Create upper and lower specification limits for 1/3 octave filter
[Lh, Ll,ff] = noctlimits(n,fc);
% Create frequency response of 1/3 octave filter with the specific cuttoff
% and sampling frequencies used here
[B,A]=noctfilt(n,fc,fs);
[H,f]=freqz(B,A,1024,fs);
% Remove zero freq. as it causes problems with log x axes
f=f(2:end);
H=H(2:end);
% Plot the FRF in magnitude (dB) and phase
figure
semilogx(f,db20(H),ff,(Lh),ff,(Ll),'LineWidth',LineWidth)
axis([fc/3 fc*3 -80 5])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Magnitude of filter char. [dB]','FontName',FontName,'FontSize',FontSize)
title('Note that the freq. axis has been zoomed in!','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
set(gca,'XTick',[200:100:1000])
% set(gca,'XTickLabel',[200,,,500,,,,,1000])
grid


%--------------------------------------------------
% Next, we create a random signal simulating a microphone signal, we then 
% A-weight this signal, and filter it using the 500 Hz third-octave filter,
% and then apply an acoustic RMS detection with a slow time constant, to
% emulate the operation of a sound level meter (SLM)
T=20;                   % Seconds long
x=randn(fs*T,1);
% First, A-weight the signal
x=timeweight(x,fs,'A');
% Then, filter it with the third-octave filter
y=filter(B,A,x);
% Then we apply an 'analog' integration filter with a slow time constant,
% i.e. with a time constant of 1 second
tau=1;                      % Time constant in s
% fci=1/(2*pi*tau);           % Cutoff freq. in Hz of integration filter
% [b,a]=butter(1,fci/(fs/2)); % Integrator filter
% y=filter(b,a,y.^2);         % Filter squared A-weighted signal
% y=y/(2*pi*fc);              % Scaled, integrated square
% y=sqrt(y);                  % Root mean square complete
y=arms(y,fs,tau);
% Now decimate this RMS(t) signal with, say, a value every half seconds
fsout=2;
n=round(fs/fsout);          % Number of samples to skip
yout=y(1:n:end);            % Decimated RMS value
t=makexaxis(yout,1/fsout);
Lp=20*log10(yout/2e-5);     % Sound pressure level in dB

% Plot this A-weighted sound level
% Note the onset of the filter during first second or so, where the RMS level
% is exponentially growing up to the approximately constant level after
% that.
figure
plot(t,Lp,'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('A-weighted SPL','FontName',FontName,'FontSize',FontSize)
title('1/3 octave at 500 Hz, Slow time constant. ','FontName',FontName,'FontSize',FontSize)