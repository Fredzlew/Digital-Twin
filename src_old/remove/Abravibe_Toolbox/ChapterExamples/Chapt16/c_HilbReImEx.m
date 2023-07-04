% c_HilbReImEx    Example of calculation of real part of FRF from imaginary
%                 and vice versa, using the Hilbert transform. It is also
%                 shown how a more causal FRF can be created
% 

% We use a measured FRF which has a non-causal impulse response, which is
% rather often seen.

% This example is similar to that used for Figure 16.6 in Brandt

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
close all
clc

%--------------------------------------------------
FontSize=11;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------
% Load an FRF from an impact test on a plate
FileName='..\data\PlexiFRF.mat';
load(FileName)
% create a new FRF H1, with imaginary part coming from Hilbert of real part
H1=frfre2im(H);
% and H2 with real part coming from imaginary...
H2=frfim2re(H);
% HH=real(H2)+j*imag(H1);
% Calculate the impulse response of the measured FRF
[h,t]=frf2ir(H,f);
% Plot results
figure;
subplot(2,2,1)
semilogy(f,abs(H),LineType{1},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Magnitude FRF [(m/s^2)/N]','FontName',FontName,'FontSize',FontSize)
grid
xlim([0 600])
title('Measured FRF','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
plot(t,h,LineType{1},'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Impulse Response [(m/s^2)s/N]','FontName',FontName,'FontSize',FontSize)
grid
axis([0 t(end) -4e5 4e5])
title('Non-causal IR','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
plot(f,real(H),LineType{1},f,real(H2),LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Real FRF [(m/s^2)/N]','FontName',FontName,'FontSize',FontSize)
grid
xlim([0 600])
title('Real of measured and from imaginary part','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
plot(f,imag(H),LineType{1},f,imag(H1),LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Imag. FRF [(m/s^2)/N]','FontName',FontName,'FontSize',FontSize)
xlim([0 600])
grid
title('Imaginary of measured and from real part','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)

% Next we calculate the IR using H1 and H2 from above, and compare
[h1,t]=frf2ir(H1,f);
[h2,t]=frf2ir(H2,f);
figure
subplot(2,1,1)
plot(t,h1,LineType{1},'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Impulse Response [(m/s^2)s/N]','FontName',FontName,'FontSize',FontSize)
grid
axis([0 t(end) -4e5 4e5])
title('IR with imaginary part of FRF replaced','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,1,2)
plot(t,h,LineType{1},'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Impulse Response [(m/s^2)s/N]','FontName',FontName,'FontSize',FontSize)
grid
axis([0 t(end) -4e5 4e5])
title('IR with real part of FRF replaced','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
