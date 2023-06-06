%[Y, df] = fftseg(y, N, dt)
%Segment a data matrix using Hanning window with 50 % overlap.
%y:     Is the data matrix, measurement channels arranged in the rows
%N:     Number of data points in each data segment
%dt:    Sampling time step
%Y:     Returned segmented data Fourier transformed
%df:    Frequency resolution
%The returned segmented data is transposed to match the preferences of the
%fft algorithm, thus in Y(r,c,s), r is the frequency index, c is the DOF
%index (1:Nc), and finally s is the segment index.

%All rights reserved, Rune Brincker, May 2012, Jul 2012.

function [Y, df] = fftseg(y, N, dt)
df = 1/(N*dt);
[nr, nc] = size(y);
ns = fix(2*nc/N);
W = [];
W1 = [0;hanning(N)];
for c=1:nr,
    W = [W, W1];
end
Y = zeros(N/2+1, nr, ns-1);
for s=1:ns-1,
    s1 = (s-1)*(N/2)+1;
    s2 = (s+1)*N/2;
    y1 = y(:,s1:s2)';
    Y2 = fft(y1.*W);
    Y(:,:,s) = Y2(1:N/2+1, :);
end


