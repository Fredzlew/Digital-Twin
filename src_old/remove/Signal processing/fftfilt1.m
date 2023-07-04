%yf = fftfilt1(y, N, dt, [optional])
%This function performs FFT filtering using a Hanning window with 50 %
%overlap. 
%y:     Signal to be filtered, measurement channels arranged in the rows
%N:     Number of data points in the segments
%dt:    Sampling time step
%yf:    Filtered signal
%The optional arguments all contain a string describing the type of 
%processing to be performed and herafter some parameters that define how. 
%The optional arguments that might appear in any order are the following:
%
%Type:  'detrend',          parameters:    mean1
%Amplify the mean in each data segment with the factor mean1
%
%Type:  'decimate',         parameters:    Ndec
%Decimate (downsample) by the factor Ndec
%
%Type:  'lowpass',          parameters:    [f1, f2]
%Lowpas filter with tapering from unity to zero from f1 to f2
%
%Type:  'highpass',         parameters:    [f1, f2]
%Highpas filter with tapering from zero to unity from f1 to f2
%
%Type:  'bandpass',         parameters:    [f0, B1, B2]
%Flat from f0-B1/2 to f0+B1/2; zero outside f0-B1/2-B2 to f0+B1/2+B2
%
%The detrending is performed by forcing the DC value of the Fourier 
%transformed data to zero. 
%All filters use a Hanning tapering in the frequency domain. If no optional
%arguments are defined the processing is neutral, except for the tapering 
%that is applied in the beginning and the end of the signal corresponding
%to the half size of the applied Hanning time window.
%
%All rights reserved, Rune Brincker, May 2012.

function yf = fftfilt1(y, N, dt, varargin)
nvar = length(varargin);
n = 1;
[Y, df] = fftseg(y, N, dt);
[nr, nc, ns] = size(Y);
H0 = ones(nr,1);
f = (0:nr-1)'*df;
om = 2*pi*f;
df = 1/(N*dt);
while (n < nvar),
    if strcmp(varargin(n),'decimate'),
        Ndec = varargin{n+1};
        nr = (nr-1)/Ndec + 1;
        Y = Y(1:nr,1:nc,1:ns)/Ndec;
        Y(nr,:,:) = real(Y(nr,:,:));
        H0 = ones(nr,1);
        f = (0:nr-1)*df;
        om = 2*pi*f;
        df = df*Ndec;
        n = n+2;
        if (n >= nvar), break; end
    end
    if strcmp(varargin(n),'detrend'),
        mean1 = varargin{n+1};
        H0(1) = mean1;
        n = n+1;
        if (n >= nvar), break; end
    end
    if strcmp(varargin(n),'lowpass'),
        f1 = varargin{n+1}(1);
        f2 = varargin{n+1}(2);
        H1 = mklpw(nr, df, f1, f2)';
        H0 = H0.*H1;
        n = n+2;
        if (n >= nvar), break; end
    end
    if strcmp(varargin(n),'highpass'),
        f1 = varargin{n+1}(1);
        f2 = varargin{n+1}(2);
        H1 = mkhpw(nr, df, f1, f2)';
        H0 = H0.*H1;
        n = n+2;
        if (n >= nvar), break; end
    end
    if strcmp(varargin(n),'bandpass'),
        f0 = varargin{n+1}(1);
        B1 = varargin{n+1}(2);
        B2 = varargin{n+1}(3);
        H1 = mkbpw(nr, df, f0, B1, B2)';
        H0 = H0.*H1;
        n = n+2;
        if (n >= nvar), break; end
    end
end
Yf = zeros(nr, nc, ns);
H = [];
for c=1:nc,
    H = [H, H0];
end
for s=1:ns,
    Yf(:,:,s) = Y(:,:,s).*H;
end
[yf, dtf] = fftsyn(Yf, df);


% W = mklpw(nf, df, f1, f2)
%Makes lowpass frequency domain window for FFT filtering using a Hanning
%tapering.
%nf:    Number of frequency lines (from DC to Nyquist)
%df:    Frequency resolution
%f1:    start Hanning tapering from unity to zero
%f2:    finish Hanning tapering
function W = mklpw(nf, df, f1, f2)
W = ones(1, nf);
f = (0:nf-1)*df;
for n=1:nf,
    if (f(n) > f1 & f(n) < f2),
        W(n) = W(n) + (cos(pi*(f(n)-f1)/(f2-f1)) - 1)/2;
    elseif (f(n) >= f2),
        W(n) = 0;
    end
end

% W = mkhpw(nf, df, f1, f2)
%Makes highpass frequency domain window for FFT filtering using a Hanning
%tapering.
%nf:    Number of frequency lines (from DC to Nyquist)
%df:    Frequency resolution
%f1:    start Hanning tapering from zero to unity
%f2:    finish Hanning tapering
function W = mkhpw(nf, df, f1, f2)
W = ones(1, nf) - mklpw(nf, df, f1, f2);

% W = mkbpw(nf, df, f0, B1, B2)
%Makes bandpass frequency domain window for FFT filtering using a Hanning
%tapering.
%nf:    Number of frequency lines (from DC to Nyquist)
%df:    Frequency resolution
%f0:    Centre frequency
%B1:    bandwidth of the flat top
%B2:    bandwidth of tapering on each side of the flap top.
function W = mkbpw(nf, df, f0, B1, B2)
W1 = mkhpw(nf, df, f0-B1/2-B2, f0-B1/2);
W2 = mklpw(nf, df, f0+B1/2, f0+B1/2+B2);
W = W1.*W2;

% W = mkbsw(nf, df, f0, B1, B2)
%Makes bandstop frequency domain window for FFT filtering using a Hanning
%tapering.
%nf:    Number of frequency lines (from DC to Nyquist)
%df:    Frequency resolution
%f0:    Centre frequency
%B1:    bandwidth of the flat bottom
%B2:    bandwidth of tapering on each side of the flap bottom.
function W = mkbsw(nf, df, f0, B1, B2)
W = ones(1, nf) - mkbpw(nf, df, f0, B1, B2);


