% R = dircor(y, N)
%Calculates the correlation function (CF) matrix by the direct method from 
%the response matrix y. The estimated correlation functions are unbiased.
%The CF matrices at different time lags are returned in the 3D array
%R(r,c,k) where r is the row entry, c is the column entry and k is the 
%discrete time lag so that the CF matrix at time lag k*dt is R(:,:,k).
%y:     Response matrix with the responses arranged in the rows.
%N:     Maximum number of time lags
%R:     Returned CF matrices

%All rights reserved. Rune Brincker, Aug 2011, July 2012.


function R = dircor(y, N)
[nr, nc] = size(y);
R = zeros(nr,nr,N/2+1);
for k=1:N/2+1,
    y1 = y(:,1:nc-k+1);
    y2 = y(:,k:nc);
    R(:,:,k) = y1*y2'/(nc-k);
end
if (nargout == 0),
    clf
    R11 = reshape(R(1,1,:),N/2+1,1);
    plot(R11)
    figure(gcf)
end
