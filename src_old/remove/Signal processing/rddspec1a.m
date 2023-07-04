% [f, Gy] = rddspec1a(y, N, dt)
%Estimates the spectral density matrix as a function of frequency based on
%random decrement (RD) functions.
%y:     Signal to be processed, measurement channels arranged in the rows    
%N:     Number of data points in the data segments (window size)
%dt:    Sampling time step
%f:     Vector containing the positive frequency axis
%Gy:    Spectral density matrices for the positive frequency axis
%The function is using the routine rddc.m for calculation of the RD
%functions.If no output arguments are used in the call, %then the singular 
%values of the spectral density matrices are plotted - or the auto spectral 
%density in case of uni-variate signal.
%Example plotting the singular values of the Heritage court data set 1:
%   load dherita1.asc
%   y = dherita1';
%   rddspec1a(y, 1024, 0.025);

%All rights reserved, Rune Brincker, May 2012.

function [f, Gyy] = rddspec1a(y, N, dt)
[nr, nc] = size(y);
df = 1/(N*dt);
nf = N/2 + 1;
f = (0:nf-1)'*df;
Gyy = zeros(nr, nr, nf);
for c=1:nr,
    [dp, dm, Np, Nm] = rddc(y, c, [-2 2], N);
    ds = dp/Np - dm/Nm;
    d = [ds(:,N/2+1:N+1), ds(:,2:N/2)];
    D = fft(d.').';
    for k=1:nf,
        Gyy(:,:,k) = Gyy(:,:,k) + conj(D(:,k))*D(:,k).';
    end
end
Gyy = Gyy/nr;
if (nargout == 0)
    if (nc == 1),
        plot(f, 10*log10(Gyy'), 'k', 'LineWidth',2)
        title('Auto spectral density')
    else,
        for k=1:nf,
            [U,S,V] = svd(Gyy(:,:,k));
            Gs(k, :) = diag(S)';
        end
        plot(f, 10*log10(Gs), 'k', 'LineWidth',2)
        title('Singular values of spectral matrix')
    end
    ylabel('dB rel. to unit')
    xlabel('Frequency [Hz]')
    figure(gcf)
    grid
end
