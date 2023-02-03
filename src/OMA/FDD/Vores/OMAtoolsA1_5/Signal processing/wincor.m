% R1 = wincor(R, W)
%Multiplies the one-sided window W onto the one-sided correlation function 
%matrix R. The correlation function matrix is a 3D array R(r,c,k) where r 
%is the row entry, c is the column entry and k is the discrete time lag so 
%that the correlation function matrix at time lag k*dt is R(:,:,k).
%If no arguments are returned then the windowed auto correlation function
%R1(1,1,:) is plotted.

%All rights reserved. Rune Brincker, May 2012.

function R = wincor(R1, W1)
[nr, nc, nk] = size(R1);
R = R1*0;
R11 = zeros(nk,1);
for r = 1:nr,
    for c = 1:nc,
        R11 = reshape(R1(r,c,:),nk,1);
        R(r,c,:) = R11.*W1;
    end
end
if (nargout == 0),
    clf
    R11 = reshape(R(1,1,:),nk,1);
    plot(R11)
    figure(gcf)
end