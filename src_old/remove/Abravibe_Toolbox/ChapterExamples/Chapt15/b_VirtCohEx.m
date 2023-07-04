% b_VirtCohEx     Calculate virtual coherences etc. for synthesized data
% 

% This example is similar to that used for Figure 15.12 in Brandt

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
% Virtual input coherence
% Define 2 SDOF systems:
fr1=100;
zr1=0.05;
fr2=2*fr1;
zr2=zr1;
m=1;
k1=(2*pi*fr1)^2*m;
k2=(2*pi*fr2)^2*m;
c1=2*zr1*sqrt(m*k1);
c2=2*zr2*sqrt(m*k1);
Type='d';

% Set up
fs=2048;
N=2*1024;
M=100;                  % Number of averages
% Add 80 dB SNR noise
NSR=1e-4;
% Correlated input signals
x3=1.5*randn(M*N,1);
x4=0.25*x3+randn(M*N,1);
% Generate outputs
y3=timefresp(x3,fs,m,c1,k1,1,1,Type);
y4=timefresp(x4,fs,m,c2,k2,1,1,Type);
yt2=y3+y4;
n2=sqrt(NSR)*std(yt2)*randn(length(yt2),1);
yt2=yt2+n2;
Gy3=apsdw(y3,fs,ahann(N));
Gy4=apsdw(y4,fs,ahann(N));
Gnn2=apsdw(n2,fs,ahann(N));
% Estimate spectra and other relationships
[Gxx2,Gyx2,Gyy2,f2]=time2xmtrx([x3 x4],yt2,fs,N);
[H2,Cm2,Cx2]=xmtrx2frf(Gxx2,Gyx2,Gyy2);
[Pca2,VCxx2]=apcax(Gxx2);

% Plot cumulated virtual coherences
figure;
subplot(1,2,1)
plot(f2,squeeze(VCxx2(:,1,1)),LineType{2},...
    f2,squeeze(VCxx2(:,1,2)),LineType{1},'LineWidth',LineWidth)
ylabel('Cum. Virt. Coh, x  _1','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 0 1.1])
grid
title('Input Signal x_1','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(1,2,2)
plot(f2,squeeze(VCxx2(:,2,1)),LineType{2},...
    f2,squeeze(VCxx2(:,2,2)),LineType{1},'LineWidth',LineWidth)
ylabel('Cum. Virt. Coh, x  _2','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 0 1.1])
title('Input Signal x_2','FontName',FontName,'FontSize',FontSize)
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

% Cont. with same data, and compute input/output virtual properties,
% similar to Fig. 15.12
[CVCyx2, VPyyx2, VGyx2] = apcaxy(Gxx2,Gyx2,Gyy2);
figure;
subplot(2,2,1)
plot(f2,CVCyx2(:,1),LineType{2},...
    f2,CVCyx2(:,2)-CVCyx2(:,1),LineType{4},...
    f2,CVCyx2(:,2),LineType{1},'LineWidth',LineWidth)
ylabel('In/Out Virt. Coh.','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 0 1.1])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
semilogy(f2,VPyyx2(:,1),LineType{2},...
    f2,VPyyx2(:,2),LineType{4},...
    f2,Gyy2-VPyyx2(:,1)-VPyyx2(:,2),LineType{3},...
    f2,Gyy2,LineType{1},'LineWidth',LineWidth)
ylabel('Virt. PSD','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 1e-18 2e-12])
grid
set(gca,'YTick',[1e-18 1e-15 1e-12])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
semilogy(f2,Gy3,LineType{2},...
    f2,Gy4,LineType{4},...
    f2,Gnn2,LineType{3},...
    f2,Gyy2,LineType{1},'LineWidth',LineWidth)
ylabel('PSD','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 1e-18 2e-12])
grid
set(gca,'YTick',[1e-18 1e-15 1e-12])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
Guu3=Gxx2(:,1,1).*abs(H2(:,1,1)).^2;
Guu4=real(Gxx2(:,2,2)).*abs(H2(:,1,2)).^2;
semilogy(f2,Guu3,LineType{2},...
    f2,Guu4,LineType{4},...
    f2,abs(Gyy2-Guu3-Guu4),LineType{3},...
    f2,Gyy2,LineType{1},'LineWidth',LineWidth)
ylabel('Coherent PSD','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
axis([0 400 1e-18 2e-12])
grid
set(gca,'YTick',[1e-18 1e-15 1e-12])
set(gca,'FontName',FontName,'FontSize',FontSize)


