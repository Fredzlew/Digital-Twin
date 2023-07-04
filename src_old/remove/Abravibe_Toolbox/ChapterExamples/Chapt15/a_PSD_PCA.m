% a_PSD_PCA     Calculate PSDs and principal components
% 

% This example is similar to that used for Figure 15.1 and 15.2 in Brandt
% (but uses somewhat different data)

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
LineType={'-','--','-.',':'};
%--------------------------------------------------

% Setup
N=2*2048;
M=100;
fs=3000;
fmin=20;
fmax=700;
NSR=1e-12;
Type='a';
FileName='..\data\PlexiSyntData';
load(FileName)
% Calculate PCAs
[Gyy,f]=acsdw(y,y,fs,N);
[Pca,CVCx]=apcax(Gyy);
for  n=1:3
    Gyy(:,n,n)=real(Gyy(:,n,n));
end
figure;
subplot(2,1,1)
semilogy(f,Gyy(:,1,1),LineType{1},...
    f,Gyy(:,2,2),LineType{1},...
    f,Gyy(:,3,3),LineType{1},'LineWidth',LineWidth);
axis([fmin fmax 1e-6 1e4])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Acceleration PSD [(m/s^2)^2/Hz]','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,1,2)
semilogy(f,Pca(:,1:3),LineType{1},'LineWidth',LineWidth);
axis([fmin fmax 1e-4 1e5])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Principal Components','FontName',FontName,'FontSize',FontSize)
grid
axis([fmin fmax 1e-6 1e4])
set(gca,'FontName',FontName,'FontSize',FontSize)
