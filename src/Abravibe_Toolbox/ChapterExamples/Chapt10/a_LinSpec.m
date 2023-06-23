% a_LinSpec    Example of computing a scaled linear spectrum

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
%
% Copyright 2011, Anders Brandt.

%
fs=1000;
N=2048;
t=(0:1/fs:(N-1)/fs)';   % Time axis
f=(0:fs/N:fs-fs/N)';	% f-axis up to sample freq, fs
f2=(0:fs/N:fs/2)';	% f-axis up to half fs (Nyquist frequency)

% Create a sine wave, amplitude 5 Volts, frequency=64 Hz
y=5*sqrt(2)*sin(2*pi*64*t);
yw=y.*hanning(N);           % Windowed time signal
Y1=abs(fft(yw));            % FFT of windowed signal
Y2=Y1/length(Y1);           % Scaling for length
Y3=winacf(hanning(N))*Y2;   % Amplitude correction
Y4=abs(Y3).^2;              % Magnitude squared
Y5=2*Y4(1:N/2+1);    % Single sided spektrum
Y6=sqrt(Y5);                % Linear spektrum

% Plot results
close
plot(t,y); title('Time data, sine, RMS=5, N=2048');
pause;
plot(t,yw); title('Windowed time data, sine, RMS=5');
pause;
plot(f,Y1); title('Magnitude of FFT of time data');
pause;
plot(f,Y2); title('FFT scaled for N');
pause;
plot(f,Y3); title('Amplitude corrected');
pause;
plot(f,Y4); title('Magnitude squared');
pause;
plot(f2,Y5); title('Single-sided autopower spectrum [V^2 RMS]');
pause;
plot(f2,Y6); title('Single-sided linear spectrum [V RMS]');
pause;

% Finally plot an overview of all steps in one plot
close
subplot(8,1,1)
plot(t,y); ylabel('Time data');
axis tight
subplot(8,1,2)
plot(t,yw); ylabel('Windowed');
axis tight
subplot(8,1,3)
plot(f,Y1); ylabel('abs(FFT)');
subplot(8,1,4)
plot(f,Y2); ylabel('FFT/N');
subplot(8,1,5)
plot(f,Y3); ylabel('Ampl. Corr.');
subplot(8,1,6)
plot(f,Y4); ylabel('abs^2');
subplot(8,1,7)
plot(f2,Y5); ylabel('Single-sided');
subplot(8,1,8)
plot(f2,Y6); ylabel('Linear');
