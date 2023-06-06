%[y, dt] = fftsyn(Y, df);
%Synthesize time signal from time segmented daya stored in the frequency
%domain (FD). Data in the FD is only stored for positive frequencies, and 
%time segmented data are assumed to be created with 50 % overlap.
%Y:     Time segmented data stored in the frequency domain.
%df:    Frequency resolution
%y:     Synthesized time signal
%dt:    Sampling time step
%Example: The following commands segments and then afterwards synthesize 
%the same signal from the Heritage Court data set.
%>>load dherita1.asc
%>>y = dherita1';
%>>N = 1024;
%>>dt = 0.025;
%>>[Y, df] = fftseg(y, N, dt);
%>>[y1, dt1] = fftsyn(Y, df);

%All rights reserved, Rune Brincker, May 2012, Jul 2012.

function [y, dt] = fftsyn(Y, df);
[nr,nc,ns] = size(Y);
N = (nr-1)*2;
dt = 1/(N*df);
Ntot = (nr-1)*(ns+1);
y = zeros(nc, Ntot);
for s=1:ns,
    s1 = (s-1)*(N/2)+1;
    s2 = (s+1)*N/2;
    ys = ifft([Y(:,:,s); flipud(conj(Y(2:N/2,:,s)))]);
    y(:,s1:s2) = y(:,s1:s2) + ys';
end
