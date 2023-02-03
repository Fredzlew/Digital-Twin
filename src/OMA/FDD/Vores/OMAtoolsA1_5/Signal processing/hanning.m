% w = hanning(N)
%Calculates a Hanning window with N-1 points, the window is symmetrical 
%around the mid point that has unit value.

%All rights reserved. Rune Brincker, Jul 2012.

function w = hanning(N)
n = N-1;
w = .5*(1 - cos(2*pi*(1:n)'/(n+1)));