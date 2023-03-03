% e_HighPassFilter      This example shows how to create a first order
%                       Butterworth highpass filter and its frequency response. 
%                       This is often a good illustration of sensor time constants, see also Chapter 7.


% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.


close all
clear

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------

%================================================================
% Settings
fs=100;                 % Sampling frequency
fc=5;                    % Cuttoff frequency of HP filter
% Create filter coefficients
fcr=fc/(fs/2);          % MATLAB/Octave uses a relative frequency axis from
                        % where 0 to 1 means 0 Hz to fs/2 Hz
[B,A]=butter(1,fcr,'high');
% Create frequency response
[H,f]=freqz(B,A,1024,fs);
% Remove zero freq. as it causes problems with log x axes
f=f(2:end);
H=H(2:end);
% Plot the FRF in magnitude (dB) and phase
subplot(2,1,1)
semilogx(f,db20(H),'LineWidth',LineWidth)
xlim([.5 fs/2])
ylabel('Magnitude of filter char. [dB]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
grid
subplot(2,1,2)
semilogx(f,angledeg(H),'LineWidth',LineWidth)
xlim([.5 fs/2])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Phase of filter char.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
grid
