% f_PSD2TimeEx  Example of generation of a time signal with known PSD

% This example is similar to that used for Figure 16.8 in Brandt

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

% Produce a wanted PSD as abs(H)^2 where H is an SDOF FRF
m=1;
k=(2*pi*100)^2*m;
c=2*0.01*sqrt(m*k);
fs=1000;
N=2048;
f=(0:fs/N:fs/2)';
H=mck2frf(f,m,c,k,1,1,'d');
Gxx=abs(H).^2;
fs=2*f(end);
% Now produce a time signal with a PSD of Gxx and with fs
L=200*1024;         % 100 blocks of 2048 samples
x=psd2time(Gxx,f,L);
% From this signal x, compute the new PSD
[Pxx,f]=apsdw(x,fs,N);
% Compute the probability density to ensure Gaussian
[PDF,Xax,GPDF]=apdf(x,33,0);
% Plot results
figure;
subplot(2,1,1)
semilogy(f,Gxx,LineType{1},f,Pxx,LineType{2},'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Displacement PSD [m^2/Hz]','FontName',FontName,'FontSize',FontSize)
legend('Target PSD','Estimated PSD')
grid
axis([0 200 1e-12 1e-7 ])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,1,2)
bar(Xax,PDF,'FaceColor',[0.8 0.8 0.8])
hold on
plot(Xax,GPDF,LineType{1},'LineWidth',LineWidth)
hold off
xlabel('Displacement [m]','FontName',FontName,'FontSize',FontSize)
ylabel('Probability Density','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

