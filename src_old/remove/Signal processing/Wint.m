% [W, f] = Wint(fs, df, f1, f2, m)
%Defines a frequency domain band-pass filter with the filter frequencies
%f1, f2. The shape of the filter is cosine tapering lifted to the exponent
%m, thus if shaping constant m=1, the tapering is pure cosine, for higher 
%values the filter becomes sharper, and for lower values of m the filter 
%becomes softer. 
%fs:       Sampling frequency
%df:       Frequency resolution
%f1, f2:   Filter frequencies
%m:        Shaping constant

%All rights reserved. Rune Brincker, Aug 2011, July 2012.

function [W, f] = Wint(fs, df, fc1, fc2, m)
fmax = fs/2;
f = 0:df:fmax;
N = length(f);
if (fc1 <= 0),
    x1 = 0;
else,
    x1 = pi*(0:df:fc1)/fc1;
end
if (fc2 >= fmax),
    x2 = 0;
else,
    x2 = pi*((fc2:df:fmax) - fc2)/(fmax - fc2);
end
f2 = fc2:df:fmax;
W1 = (1 - cos(x1))/2;
W2 = (1 + cos(x2))/2;
N1 = length(W1);
N2 = length(W2);
W0 = ones(1, N - N1 - N2);
W = [W1, W0, W2];
W = W.^m;