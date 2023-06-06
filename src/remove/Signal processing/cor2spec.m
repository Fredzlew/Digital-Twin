% [f, G] = cor2spec(R, dt)
%Calculates the spectral density (SD) matrix G from the corelation function 
%(CF)matrix R. Both G and R are 3D arrays G(r,c,k), R(r,c,k) where r is the 
%row entry, c is the column entry and where k for the SD matrix is the 
%discrete frequency so that the SD matrix at frequency k*df is G(:,:,k) and 
%for the CF matrix it is the discret time lag so that the CF matrix at time 
%lag k*dt is R(:,:,k). The parameter dt is the time step and f a vector
%holding the frequencies in Hz.

%All rights reserved. Rune Brincker, May 2012.

function [f, G1] = cor2spec(R1, dt)
[nr, nc, nf] = size(R1);
N = (nf-1)*2;
T = N*dt;
df = 1/T;
G1 = R1*0;
R11 = zeros(nf,1);
for r=1:nr,
    for c=1:nc,
        R11 = reshape(R1(r,c,:),nf,1);
        R = [R11; flipud(R11(2:nf-1))];
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
