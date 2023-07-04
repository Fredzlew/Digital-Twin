% W = flattrian(N, a)
%Defines and returns the flat triangular window. The total number of points
%is N-1 so that the window is symmetric around the mid point. The paramter
%a defines where to change from the flat part (constant window) to the
%linear part that decreases to zero at the boundaries. the window is
%returned as a column vector. Parameter a = 0 defines a triangular window,
%and a = 1 defines a flat (boxcar) window.

%All rights reserved. Rune Brincker, May 2012.

function W = flattrian(N, a)
W1 = zeros(N/2,1);
k1 = a*N/2;
for k=1:N/2,
    if (k <= k1),
        W1(k) = 1;
    else,
        W1(k) = (N/2 - k + 1)/(N/2 - k1 + 1);
    end
end
W = [flipud(W1(2:N/2)); W1];