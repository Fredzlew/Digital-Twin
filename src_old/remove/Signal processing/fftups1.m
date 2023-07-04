% y1 = fftups1(y, Nups)
% This function upsample the data to sampling rate that is Nups times 
% higher than the original sampling rate.
% y:    Signal matrix, signals arranged in the rows
% Nups: Number of times to upsample
% Example: ex4_1

% All rights reserved, Rune Brincker, Aug 2011, July 2012.

function y1 = fftups1(y, Nups)
y = y';
[N, Nc] = size(y);
Z = zeros(N*(Nups-1), Nc);
Y = fft(y);                           %Step 1, take the Fourier transform
Y1 = [Y(1:N/2+1,:); Z; Y(N/2+2:N);];  %Step 2, pad zeroes in the frequency domain
y1 = Nups*real(ifft(Y1))';             %Step 3, Transform back
t = 0:(N-1);
t1 = (0:N*Nups-1)/Nups;
if (nargout == 0),
    clf
    plot(t1, y1(:,1), 'r*')
    hold on
    plot(t, y(:,1), 'o')
    figure(gcf)
end

