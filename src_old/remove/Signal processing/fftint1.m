% y = fftint1(x, fs, f1, f2, m, n)
% 
% This function integrates/differentiates a signal matrix using one-shot 
% FFT-filtering.
% No windowing is used on the data matrix x. 
% The function uses a band-pass frequency domain filter decribed by the 
% two filter frequencies f1, f2. 
%
% x:     Data matrix input
% fs:    Sampling frequency
% f1:	 Low cut-off frequency
% f2:	 High cut-off frequency
% m:	 Filter shape parameter (win = han.*m)
% N:	 Integration parameter, n=0 does nothing, N=1 integrates one time, 
%           n=-1 differentiates one time 
% y:      Integrated data matrix

% All rights reserved, Rune Brincker, Nov 2010, Aug 2011, July 2012.

function y = fftint1(x, fs, f1, f2, m, n)
Y = x.'*0;
[N, Nc] = size(Y);
dt = 1/fs;
T = N*dt;
df = 1/T;
Nf = N/2+1;
[W, f] = Wint(fs, df, f1, f2, m);               %Step 1, define filtering window
H1 = zeros(Nf, 1);
H1(2:Nf-1) = W(2:Nf-1).*((i*2*pi*f(2:Nf-1)).^(-n));   %Step 2, define one sided FRF
H = [H1; conj(flipud(H1(2:Nf-1)))];             %Step 2, define double sided FRF
X = fft(x.');                                   %Step 3, Calculate FFT to get X
for c = 1:Nc,
    Y(:,c) = X(:, c).*H;                        %Step 4, Multiply X and H
end
y = ifft(Y).';                                  %Step 5, Take inverse FFT






