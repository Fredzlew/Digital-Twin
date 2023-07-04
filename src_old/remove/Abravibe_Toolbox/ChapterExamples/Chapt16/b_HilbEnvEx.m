% b_HilbEnvEx     Example of calculation of envelope using Hilbert transform
% 

% This example uses narrowband noise and two different time delays to make
% a cross-correlation signal with two distinct peaks, one at each time
% delay. The peaks are found easier in the envelope.

% This example is similar to that used for Figure 16.4 and 16.5 in Brandt

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
% Define parameters
fs=2000;
N=2*2048;
% Define a SDOF systems
f1=200;
z1=0.01;
m=1;
k=(2*pi*f1)^2*m;
c=2*z1*sqrt(m*k);
% Define two different propagation delays
T1=0.1;
T2=0.15;
% Define random input noise
x=randn(500*N,1);
% Make the noise narrow band by passing it as a force to the SDOF system,
% and use the displacement output
y=timefresp(x,fs,m,c,k,1,1,'d');
% Delay y by T1 and T2 respectively, to produce the two signals y1 and y2
y1=[zeros(T1*fs,1);y];
y2=[zeros(T2*fs,1);y];
% Set the smallest length signal and cut remaining signals (axcorr requires
% all signals to be equal length)
L=length(y);
y1=y1(1:L);
y2=y2(1:L);
% Sum the two delayed signals and compute cross-correlation between the
% undelayed signal and the sum signal
ys=y1+y2;
[Ryx,tau]=axcorr(y,ys,fs,N);
% Calculate the Hilbert envelope
z=hilbert(Ryx);
e=abs(z);
% Plot results
figure;
plot(1000*tau,1e6*Ryx,LineType{4},1000*tau,1e6*e,LineType{1},'LineWidth',LineWidth);
xlabel('Time lag [ms]','FontName',FontName,'FontSize',FontSize)
ylabel('R_{yx} (\tau) and envelope','FontName',FontName,'FontSize',FontSize)
grid
xlim([0 200])
set(gca,'FontName',FontName,'FontSize',FontSize)
