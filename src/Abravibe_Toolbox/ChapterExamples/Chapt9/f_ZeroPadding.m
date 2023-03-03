% f_ZeroPadding     Test of zero padding using two closely spaced sines
%                   Does zero padding help to increase the frequency resolution?

clear
close all
clc

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% Check fft of a zero padded cosine
fs=1024;
N=1024;
Z=16;                            % Zero padding to Z*N
w=ahann(N);
t=makexaxis((1:N)',1/fs);
y1=cos(16.5*2*pi*t)+cos(17.75*2*pi*t);

% Compare linear spectra with Hanning window of both
Y1=sqrt(2)*2*abs(fft(w.*y1,N)/N);       % DFT without zero padding
Y1=Y1(1:N/2+1);
f1=(0:fs/N:fs/2)';
Y2=sqrt(2)*2*Z*abs(fft(w.*y1,Z*N)/Z/N); % DFT with strong zero padding
Y2=Y2(1:Z*N/2+1);
f2=(0:fs/Z/N:fs/2)';
figure
plot(f1,Y1,'k-',f2,Y2,'k-.')
title('Linear spectra of cosines, exactly on freq. line')
legend('No zero padd','16x Zero Padded')
axis([0 30 0 1])

