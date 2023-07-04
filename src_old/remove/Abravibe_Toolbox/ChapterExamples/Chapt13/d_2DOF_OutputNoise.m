% d_2DOF_OutputNoise    Example of SISO FRF estimation with pure, burst, and 
%                       pseudo random input noise, and output extraneous noise. 
%                       This example is similar to the 
%                       example in Figure 13.19 in Brandt. 
%

% This example uses output contaminating noise.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% 2012-01-15 Fixed bug, erroneously calling for synchronous averaging, and
%            made running better in Octave

clear 
close all
clc

%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-','--',':','o','+'};


fs=2048;
fmin=10;
fmax=500;
Npure=16384;
Nburst=4*1024;
BurstLength=25;
NPseudo=4*1024;
NoBlocks=100;
SN=1e-3;
OutputType='d';         % Displacement. DO NOT change, as the true FRF will not be correct
f=(0:fs/(64*1024):fs/2)';
fidx=find(f>10 & f < 500);
% Set up 2DOF model
fr=[100 200];
zr=0.01*[1 1];
p=-2*pi*zr.*fr+j*2*pi*fr.*sqrt(1-zr.^2);
V=[1 1;1 -1];
Htrue=modal2frf(f,p,V,1,1,OutputType);
Htrue=Htrue(fidx);
ftrue=f(fidx);
% Pure random
x=randn(NoBlocks*Npure,1);
y=timefresp(x,fs,p,V,1,1,OutputType);
n=std(y)*SN*randn(length(y),1);
y=y+n;
[Gxx,Gyx,Gyy,fpure,NBlocks1]=time2xmtrx(x,y,fs,hsinew(Npure),67);
[Hpure,Cpure]=xmtrx2frf(Gxx,Gyx,Gyy);
fpureidx=find(fpure > fmin & fpure < fmax);
Hpure=Hpure(fpureidx);
Cpure=Cpure(fpureidx);
fpure=fpure(fpureidx);
% Burst random
x=abrand(Nburst,NoBlocks,BurstLength(1));
y=timefresp(x,fs,p,V,1,1,OutputType);
y=y+n(1:length(y));
[Gxx,Gyx,Gyy,fburst,NBlocks2]=time2xmtrx(x,y,fs,boxcar(Nburst),0);
[Hburst,Cburst]=xmtrx2frf(Gxx,Gyx,Gyy);
fburstidx=find(fburst > fmin & fburst < fmax);
Hburst=Hburst(fburstidx);
Cburst=Cburst(fburstidx);
fburst=fburst(fburstidx);
%Pseudo random
NPBlocks=NoBlocks;
x=aprand(NPseudo);
xt=x;
for n=1:NPBlocks
    x=[x;xt];
end
y=timefresp(x,fs,p,V,1,1,OutputType);
y=y+std(y)*SN*randn(length(y),1);
% Throw away transient part
x=x(11*NPseudo+1:end);
y=y(11*NPseudo+1:end);
[Gxx,Gyx,Gyy,f3,NBlocks3]=time2xmtrx(x,y,fs,boxcar(NPseudo),0);
[H3,C3]=xmtrx2frf(Gxx,Gyx,Gyy);
f3idx=find(f3>10 & f3 < 500);
Hpseudo=H3(f3idx);
fpseudo=f3(f3idx);
Cpseudo=C3(f3idx);
% Plot results
figure;
subplot(2,2,1)
semilogy(fpure,abs(Hpure),fburst,abs(Hburst),...
    fpseudo,abs(Hpseudo),ftrue,abs(Htrue),'LineWidth',LineWidth);
% axis([99 101 1e-4 1.3e-4])
xlim([fmin fmax])
ylabel('Magnitude FRF','FontName',FontName,'FontSize',FontSize)
title('All FRFs','FontName',FontName,'FontSize',FontSize)
legend('Pure','Burst','Pseudo','True')
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
plot(fpure,Cpure,fburst,Cburst,...
    fpseudo,Cpseudo,'LineWidth',LineWidth);
axis([fmin fmax 0.8 1.01])
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Coherence','FontName',FontName,'FontSize',FontSize)
title('All coherences','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
semilogy(fpure,abs(Hpure),LineType{1},fburst,abs(Hburst),LineType{4},...
    fpseudo,abs(Hpseudo),LineType{5},ftrue,abs(Htrue),LineType{2},'LineWidth',LineWidth);
title('All FRFs','FontName',FontName,'FontSize',FontSize)
legend('Pure','Burst','Pseudo','True')
axis([98 102 .1 .17])
subplot(2,2,4)
idxpure=find(fpure >= 155.5 & fpure <= 160.5);
idxburst=find(fburst >= 155.5 & fburst <= 160.5);
idxpseudo=find(fpseudo >= 155.5 & fpseudo <= 160.5);
idxtrue=find(ftrue >= 155.5 & ftrue <= 160.5);
fpure=fpure(idxpure);
Hpure=Hpure(idxpure);
fburst=fburst(idxburst);
Hburst=Hburst(idxburst);
fpseudo=fpseudo(idxpseudo);
Hpseudo=Hpseudo(idxpseudo);
ftrue=ftrue(idxtrue);
Htrue=Htrue(idxtrue);
semilogy(fpure,abs(Hpure),LineType{1},fburst,abs(Hburst),LineType{4},...
    fpseudo,abs(Hpseudo),LineType{5},ftrue,abs(Htrue),LineType{2},'LineWidth',LineWidth);
% axis([155.5 160.5 2e-4 4e-4])
% xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
% title('All FRFs','FontName',FontName,'FontSize',FontSize)

