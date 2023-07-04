% pltspec(G, dt, ns, [var])
%This function plots spectral density functions. If a 1D array is supplyed
%in the the input variable G, then this spectral density is plotted; if a
%3D array is supplyed in G, then the singular values of the spectral
%density matrix is plotted. 
%G:     Spectral density to be plotted, either 1D or 3D data
%dt:    Sampling time step
%ns:    Number of singular values to be plotted
%var:   Optional variables to specify colors etc. like in Matlabs plot 
%       function.

% All rights reserved, Rune Brincker, Jul 2012.


function pltspec(G, dt, ns, varargin)
[nr, nc, nf] = size(G);
if (nr == 1 & nc == 1),
    G = reshape(G(1,1,:), nf, 1);
    plth(G, dt, varargin{:} )
else,
    for k=1:nf,
        [U,S,V] = svd(G(:,:,k));
        Gs(k, :) = diag(S(1:ns,1:ns))';
    end
   
    df = 1/((nf-1)*2*dt);
    f = (0:nf-1)'*df;
    plot(f, 10*log10(Gs), varargin{:})
    ylabel('dB rel. to unit')
    xlabel('Frequency [Hz]')
%     title('Singular values of spectral matrix')
    grid
end

