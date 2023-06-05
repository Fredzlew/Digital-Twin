% c_SISOBurstRandEx Example of SISO FRF estimation with burst random input noise
%                   This example is similar to the example in Figure 13.18 
%                   in Brandt. It computes the FRF using two burst random
%                   percentages and compares with the true FRF
%
% Please see the discussion about accuracy in forced response measurements
% in example i_3DOF_MKz in the folder for chapter 6

% This example uses no input or output contaminating noise.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear 
close all
clc

%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

% Set up simulation parameters (three different blocksizes)
% N=1024*[2 4 16];
NoBlocks=50;
fs=2048;
fmin=10;                
fmax=500;
%==================================================================

% Set-up an SDOF mechanical system
OutputType='d';         % Displacement. DO NOT change, as the true FRF will not be correct
fr=100;                 % Undamped natural frequency
zr=.01;                 % Relative damping
m=1;                    % Mass = 1 kg
k=(2*pi*fr)^2*m;        % Stiffness to achieve fr
c=zr*2*sqrt(m*k);       % Damping to achieve zr

% Compute modal parameters, poles and "mode shape"
[p,V]=mck2modal(m,c,k);


% % Synthesize true accelerance FRF H_RespDof,RefDof using the following
% lines if you do not care about 'absolute' accuracy. It works in most
% cases as the timefresp is very accurate. The code in the plotting section
% below is prepared if you want to check what the results would be with the
% following lines instead of the default lines. Then uncomment the next
% four lines, and the two lines for plot c)
ftrue=(0:fs/(64*1024):fs/2)';
fidx=find(ftrue>fmin & ftrue < fmax);
ftrue=ftrue(fidx);
Htrue=modal2frf(ftrue,p,V,1,1,OutputType);

% Compute the response signal for a Gaussian input force
% To be able to show how exact the estimates are compared to the true
% values, we let timefresp respond with the exact FRF used by the digital
% filters, and not the FRF of the modal model. In most cases you can
% actually use the approximate FRF as showed in the outcommented lines
% above the present lines.
NBurst=2*1024;
BurstLength=[25 75];
for n=1:length(BurstLength)
    x=abrand(NBurst,NoBlocks,BurstLength(n));
    y=timefresp(x,fs,p,V,1,1,OutputType);
    [Gxx,Gyx,Gyy,f2,NBlocks2]=time2xmtrx(x,y,fs,boxcar(NBurst),0);
    [H2,C2]=xmtrx2frf(Gxx,Gyx,Gyy);
    f2idx=find(f2>10 & f2 < 600);
    Hburst{n}=H2(f2idx);
    fburst{n}=f2(f2idx);
    Cburst{n}=C2(f2idx);
end
hf=figure;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 12]);
set(gcf, 'PaperSize', [12 12])
subplot(2,2,1)
semilogy(fburst{1},abs(Hburst{1}),LineType{4},fburst{2},abs(Hburst{2}),LineType{3},...
    ftrue,abs(Htrue),LineType{1},'LineWidth',LineWidth);
axis tight
ylabel('Magnitude FRF','FontName',FontName,'FontSize',FontSize)
title('a)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
legend('25% burst','75% burst','True FRF')
subplot(2,2,3)
plot(fburst{1},Cburst{1},LineType{4},fburst{2},Cburst{2},LineType{3},...
    'LineWidth',LineWidth);
axis([fmin fmax 0 1.1])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Coherence','FontName',FontName,'FontSize',FontSize)
title('b)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
semilogy(fburst{1},abs(Hburst{1}),'ko',fburst{2},abs(Hburst{2}),'k+',... %LineType{3},...
    ftrue,abs(Htrue),LineType{1},'LineWidth',LineWidth);
axis([95 105 2e-5 2e-4])
legend('25% burst','75% burst','True FRF')
set(gca,'YTick',[])
title('c)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
plot(fburst{1},Cburst{1},LineType{4},fburst{2},Cburst{2},LineType{3},...
    'LineWidth',LineWidth);
axis([95 105 0.5 1.1])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
title('d)','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
