% W = exponential(N, eps, dt)
%Defines and returns the exponential window. The total number of points
%is N-1 so that the window is symmetric around the mid point. The paramter
%eps defines defines the value of the window at the boundaries (window 
%decays from unity to eps). the window is returned as a column vector.

%All rights reserved. Rune Brincker, May 2012.

function W = exponential(N, eps, dt)
T = N*dt;
a = (2/T)*log(eps);
W1 = zeros(N/2,1);
t = (0:N/2-1)'*dt;
W1 = exp(a*t);
W = [flipud(W1(2:N/2)); W1];
