% [dp, dm, Ntp, Ntm] = rddc(y, c, e, N)
%Calculates the random decrement (RD) signatures for a level trig condition
%where the trigering levels are defined by the paramter e the triggering 
%channel is defined by the channel number c and the maximum time lag is
%defined by N.
%y:     Signal matrix where the signals are arranged in the rows.
%c:     Channel number to apply the level trig condition
%e:     Dimensionless trig level so that the actual trig level is given by 
%       a = std(y(c,:))*e;
%N:     Maximum time lag in RD functions
%dp:    One sided RD function matrix with positive slope
%dm:    One sided RD function matrix with negative slope
%Ntp:   Number of triggering points with positive slope
%Ntm:   Number of triggering points with negative slope
%The RD function matrices are arranged with RD functions in the rows. The
%function can also be used to perform the similar RD estimates for a level
%band, in this case the the paramter e is a vector with two elements 
%[e1,e2] specifying the band interval. 

%All rights reserved. Rune Brincker, May 2012.

function [dp, dm, Ntp, Ntm] = rddc(y, c, e, N)
[nr, nt] = size(y);
a = std(y(c,:))*e;
dp = zeros(nr, N+1);
dm = zeros(nr, N+1);
na = length(a);
Nmin = N/2 + 2;
Nmax = nt - N/2;
Ntp = 0;
Ntm = 0;
if (na == 1),
    for k = Nmin:Nmax,
        if (y(c, k-1) < a & y(c, k) > a),
            dp = dp + y(:,k-N/2:k+N/2) + y(:,k-N/2-1:k+N/2-1);
            Ntp = Ntp + 1;
        elseif (y(c, k-1) > a & y(c, k) < a),
            dm = dm + y(:,k-N/2:k+N/2) + y(:,k-N/2-1:k+N/2-1);
            Ntm = Ntm + 1;
        end
    end
end
if (na == 2),
    for k = Nmin:Nmax,
        if (a(1) < y(c, k) & y(c, k) < a(2))
            if (y(c, k) > y(c, k-1)),
                dp = dp + y(:,k-N/2:k+N/2) + y(:,k-N/2-1:k+N/2-1);
                Ntp = Ntp + 1;
            else,
                dm = dm + y(:,k-N/2:k+N/2) + y(:,k-N/2-1:k+N/2-1);
                Ntm = Ntm + 1;
            end
        end
    end
end

          
            