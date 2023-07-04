% [f, Gy] = halfspecw(y, N, dt)
%Estimates the spectral density matrix as a function of frequency based on
%the direct calculated correlation function matrix and hereafter forcing
%the negative part of all correlation functions to zero. The half outermost
%part of the correlation functions are forced to zero using a flat
%triangular window.
%y:     Signal to be processed, measurement channels arranged in the rows    
%N:     Number of data points in the data segments (window size)
%dt:    Sampling time step
%f:     Vector containing the positive frequency axis
%Gy:    Spectral density matrices for the positive frequency axis
%The function is using the routine dircor.m for calculation of the 
%correlation function matrix. If no output arguments are used in the call, 
%then the singular values of the spectral density matrices are plotted - 
%or the auto spectral density in case of uni-variate signal.
%Example plotting the singular values of the Heritage court data set 1:
%   load dherita1.asc
%   y = dherita1';
%   halfspecw(y, 1024, 0.025);

%All rights reserved, Rune Brincker, May 2012.

function [f, G1] = halfspecw(y, N, dt)
R1 = dircor(y, N);
[nr, nc, nf] = size(R1);
N = (nf-1)*2;
W = flattrian(N, 0.5);
W1 = zeros(N/2+1,1);
W1(1:N/2) = W(N/2:N-1);
R2 = wincor(R1, W1);
T = N*dt;
df = 1/T;
G1 = R1*0;
R11 = zeros(nf,1);
for r=1:nr,
    for c=1:nc,
        R11 = reshape(R2(r,c,:),nf,1);
        R = [R11; flipud(R11(2:nf-1))*0];
        G = fft(R);
        G1(r,c,:) = G(1:nf);
    end
end
f = (0:nf-1)*df;
Gs = zeros(nf, nc);
if (nargout == 0)
    if (nc == 1),
        plot(f, 10*log10(G1'), 'k', 'LineWidth',2)
        title('Auto spectral density')
    else,
        for k=1:nf,
            [U,S,V] = svd(G1(:,:,k));
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