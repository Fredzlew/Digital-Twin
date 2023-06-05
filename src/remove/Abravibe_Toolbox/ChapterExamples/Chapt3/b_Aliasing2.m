% b_Aliasing2.          This example only works in MATLAB!
%
% A sine wave with frequency stepped up passed the Nyquist frequency
% is sampled at fs Hz. A Spectrum is plotted simultaneously with outputting
% the sine to the sound card of the pc. You hear the frequency going
% up continuously, while you will see the sampled signal spectrum decrease
% in frequency once it passes half the sampling frequency.


% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.


% Initialize
clear
close all

% Parameters
fs=2000;         % Sampling frequency for simulation
fs2=44100;       % fs for audio out; not all audio boards support 'odd' sampling frequencies
T=1;             % Duration of each tone

fvect=100*[5:15]; % Sine freqs from 100 to 2000 Hz

t=(0:1/fs:T)';     % Time axis for T seconds of data
t2=(0:1/fs2:T)';    % Audio sampling time axis

warning off        % Never mind this.
figure
ys=sin(2*pi*fvect(1)*t);    % Sampled sine
[Y,f]=alinspec(ys,fs,ahann(length(ys)),1,0);
plot(f,Y)
axis([0 1000 0 1])
xlabel('Frequency [Hz]')
ylabel('Spectrum')
title('Sampling freq. = 2000 Hz, half = 1000 Hz')
fprintf('Press <RETURN> to start...\n')
pause
for n=1:length(fvect)
%     subplot(1,2,1)
    ys=sin(2*pi*fvect(n)*t);    % Sampled sine
    yout=sin(2*pi*fvect(n)*t2);
%     plot(t,ys)
%     subplot(1,2,2)
    [Y,f]=alinspec(ys,fs,ahann(length(ys)),1,0);
    plot(f,Y)
    axis([0 1000 0 1])
    xlabel('Frequency [Hz]')
    ylabel('Spectrum')
    title(['fs = 2000 Hz, Nyquist = 1000 Hz, Inst. freq. = ' num2str(fvect(n))])
    soundsc(yout,fs2)
    pause(1)
end
warning on      % Never mind this either.