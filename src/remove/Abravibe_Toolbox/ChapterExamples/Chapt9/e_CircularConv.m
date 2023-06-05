% e_CircularConv    Illustrate circular convolution and zero padding
%                   Using a test sequence x=[1 1 1 1 1]

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
close all

% Test sequence(s)
x=[1 1 1 1].';
y=[1 2 3 0].';
figure
subplot(2,1,1)
stem(1:length(x),x,'filled')
ylabel('x')
title('Test sequences')
subplot(2,1,2)
stem(1:length(y),y,'filled')
ylabel('y')

% Linear convolution
zc=conv(x,y);
figure
stem(1:length(zc),zc,'filled')
title('Linear (true) convolution of x by y')

% Circular convolution by direct DFT
Pxy=fft(x).*fft(y);             % FFT product
plotreim(1:length(x),Pxy,'k','filled')
subplot(2,1,1)
title('Product X * Y')

zcc=ifft(Pxy);
plotreim(1:length(zcc),zcc,'k','stem')
subplot(2,1,1)
title('Real and imag part of circular convolution')
% This is apparently NOT like the linear correlation

% Circular convolution with zero padding (=linear convolution)
xz=[x;zeros(length(x),1)];
yz=[y;zeros(length(xz)-length(y),1)];
Pxyzz=fft(xz).*fft(yz);
plotreim(1:length(Pxyzz),Pxyzz,'k','stem')
title('Product of zero padded spectra')
zzpcc=ifft(Pxyzz);
plotreim(1:length(zzpcc),zzpcc,'k','stem')
title('Circular conv. with zero padding length(x)')


