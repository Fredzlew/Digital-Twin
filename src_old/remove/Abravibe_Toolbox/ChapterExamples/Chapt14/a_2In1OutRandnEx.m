% a_2In1OutRandnEx   Simulate 2-in-1-out system excited by random noise
% 
% This example includes output contaminating noise, but it is set to 0 by
% default. Change NSRout below to, e.g. 1e-4 to add some noise.

% This example is similar to that used for Figure 14.11 in Brandt although
% it uses pure random and not burst random as in Figure 14.11.

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

%================================================================
% 2In-1Out pure random, fixed noise level on output only
 
% Setup
fs=2048;
fmin=10;
fmax=500;
N=16*1024;
NoBlocks=50;
NSRout=1e-4;
f=(0:fs/(64*1024):fs/2)';
fidx=find(f>fmin & f < fmax);

% Set up 2DOF model in modal parameters
OutputType='d';
fr=[100 200];
zr=0.01*[1 1];
p=-2*pi*zr.*fr+j*2*pi*fr.*sqrt(1-zr.^2);
V=[1 1;1 -1];
Htrue=modal2frf(f,p,V,[1:2],1,OutputType);
Htrue=squeeze(Htrue);
Htrue=Htrue(fidx,:);
ftrue=f(fidx);
%------------------------
% Create random input signals
x1=randn(N*NoBlocks,1);
x2=randn(N*NoBlocks,1);
% Create output signal without noise
y=timefresp(x1,fs,p,V,1,1,OutputType)+timefresp(x2,fs,p,V,2,1,OutputType);
% Add output noise
n=sqrt(NSRout)*std(y)*randn(length(y),1);
y=y+n;
% Estimate MIMO system
[Gxx,Gyx,Gyy,fb,NBb]=time2xmtrx([x1 x2],y,fs,ahann(N));
[Hb,Cmb,Cxb]=xmtrx2frf(Gxx,Gyx,Gyy);
Hb=squeeze(Hb);                         % Change from Nf-by-1-by-2 to Nf-by-2
fbidx=find(fb > fmin & fb < fmax);
Hb=Hb(fbidx,:);
fb=fb(fbidx);
Cmb=Cmb(fbidx);
Cxb=Cxb(fbidx,:,:);

% Plot results
figure
subplot(3,2,1)
semilogy(fb,abs(Hb(:,1)),LineType{1},'LineWidth',LineWidth);
xlim([fmin fmax])
ylabel('|H_{11}|','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,2,2)
semilogy(fb,abs(Hb(:,1)),'ko',ftrue,abs(Htrue(:,1)),LineType{1},...
    'LineWidth',LineWidth);
axis([98.2 101.8 .08 .17])
ylabel('|H_{11}|','FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,2,3)
semilogy(fb,abs(Hb(:,2)),LineType{1},'LineWidth',LineWidth);
xlim([fmin fmax])
ylabel('|H_{12}|','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,2,4)
semilogy(fb,abs(Hb(:,2)),'ko',ftrue,abs(Htrue(:,2)),LineType{1},...
    'LineWidth',LineWidth);
axis([98.2 101.8 .08 .17])
ylabel('|H_{12}|','FontName',FontName,'FontSize',FontSize)
set(gca,'YTick',[])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,2,5)
plot(fb,Cmb,LineType{1},'LineWidth',LineWidth);
xlim([fmin fmax])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('\gamma^2_{y:x}','FontName',FontName,'FontSize',FontSize)
ylim([0 1.05])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(3,2,6)
plot(fb,Cxb(:,1,2),LineType{1},'LineWidth',LineWidth);
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
xlim([fmin fmax])
ylabel('\gamma^2_{12}','FontName',FontName,'FontSize',FontSize)
ylim([0 1.05])
set(gca,'FontName',FontName,'FontSize',FontSize)

