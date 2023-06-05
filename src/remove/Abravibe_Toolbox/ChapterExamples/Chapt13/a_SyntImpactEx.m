% a_SyntImpactEx   Impact excitation example using synthesized data

clear
clc
close all

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};

%--------------------------------------------------
% Create suitable pulse for plate impact
% To make the example realistic, we show a common effect of reduced
% sampling frequency: When the sampling frequency used to measure the
% force impulse is below twice the bandwidth of the pulse, 'ringing' will
% occur in the sampled force pulse. This is quite normal and does not pose
% a problem from a measurement point of view.
fs=2000;					% Pulse sampling frequency
fsn=500;                    % Used sampling frequency
N=2*2048;
NN=fs/fsn*N;
p1=makepulse(NN,fs,.005,'Halfsine');
p1=[p1(end-200:end);p1];
p1=p1(1:NN);
% Resample down to fsn
p1=resample(p1,fsn,fs);
%===========================
% Add a little noise
p1=p1+0.001*max(p1)*randn(N,1);
%===========================
fs=fsn;
w=aforcew(N,2,'Gaussian');
f1=p1.*w;
f2=p1; %f2(150:end)=0;
t=makexaxis(p1,1/fs);
[T1,f]=atranspec(f1,fs);
[T2,f]=atranspec(f2,fs);
% Plot in 4 quadrants, impact, force win, spectrum with/without win,
% similar to fig. 13.10
figure
subplot(2,2,1)
plot(t,f2,LineType{1},'LineWidth',LineWidth)
axis([0 5 -.01 .4])
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Force signal [N]','FontName',FontName,'FontSize',FontSize)
grid
title('Force; zoom in to see ringing','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
plot(t,w,LineType{1},'LineWidth',LineWidth)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Force window','FontName',FontName,'FontSize',FontSize)
grid
axis([0 5 -.01 1.1])
title('Force window','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,2)
semilogy(f,T2,LineType{1},'LineWidth',LineWidth)
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Transient Spectrum [Ns]','FontName',FontName,'FontSize',FontSize)
title('Unwindowed spectrum','FontName',FontName,'FontSize',FontSize)
xlim([0 200])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
semilogy(f,T1,LineType{1},'LineWidth',LineWidth)
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Transient Spectrum [Ns]','FontName',FontName,'FontSize',FontSize)
title('Windowed spectrum','FontName',FontName,'FontSize',FontSize)
xlim([0 200])
set(gca,'FontName',FontName,'FontSize',FontSize)


%==========================================================================
% Next we let the force pulse excite a structure simulated by modal
% properties from a slalom ski (see example b_ImpactMeasEx in Chapt 13)
% And then add some noise. An exponential window is added and spectra of
% unwindowed and windowed responses are plotted
amin=1e-9;amax=1e-5;
RefDof=14;
RespDof=1;
% Load modal model and put mode shapes and poles into V and p, respectively
% Model has only Z dofs
load('..\data\SkiModes')
% Synthesize true accelerance FRF H_RespDof,RefDof
f=(0:fs/8/N:fs/2)';
fidx=find(f>10 & f < 600);
Htrue=modal2frf(f,p,V,RefDof,RespDof,'d');
% Now compute the response signal for the pulse porg
y=timefresp(p1,fs,p,V,RefDof,RespDof,'d');
%===========================
% Add some output noise to the response signal
y=y+max(abs(y))/100000*randn(length(y),1);
%===========================
% Windowed version
w=aexpw(length(y),.001);
yw=y.*w;
t=makexaxis(y,1/fs);        % For plotting
% Calculate spectra of responses
[Ty,f]=atranspec(y,fs);
[Tyw,f]=atranspec(yw,fs);
figure
subplot(2,2,1)
plot(t,y,LineType{1},'LineWidth',LineWidth)
title('Driving point acceleration, no exp. window','FontName',FontName,'FontSize',FontSize)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Displacement [m]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
% axis([0 t(end) -2.5 2.5])
grid
subplot(2,2,3)
plot(t,yw,LineType{1},t,w*max(y),LineType{4},'LineWidth',LineWidth)
title('Driving point acceleration, exp. window w(end)=.001','FontName',FontName,'FontSize',FontSize)
xlabel('Time [s]','FontName',FontName,'FontSize',FontSize)
ylabel('Displacement [m]','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
% axis([0 t(end) -2.5 2.5])
grid
subplot(2,2,2)
semilogy(f,Ty,LineType{1},'LineWidth',LineWidth)
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Transient Spectrum [ms]','FontName',FontName,'FontSize',FontSize)
title('Unwindowed spectrum','FontName',FontName,'FontSize',FontSize)
axis([0 200 amin amax])
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
semilogy(f,Tyw,LineType{1},'LineWidth',LineWidth)
grid
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Transient Spectrum [ms]','FontName',FontName,'FontSize',FontSize)
title('Windowed spectrum','FontName',FontName,'FontSize',FontSize)
axis([0 200 amin amax])
set(gca,'FontName',FontName,'FontSize',FontSize)

%=========================================================================
% Now simulate a measurement with 5 averages and equal pulses (only the
% extraneous noise is changing between the averages)
fmin=0;fmax=200;
amin=1e-7;amax=1e-2;
% Set the analysis parameters here:
%-----------------------------------
NumberAverages=5;
NSRForce=1e-4;      % Linear (amplitude) factor
NSRResponse=2e-4;
ExpWinEnd=1;     % Note! In percent
%-----------------------------------
% Simulate the averaging process by using the pulse p1 plus random noise,
% then add a window to each block, and totally produce 5 blocks in pout and
% yout for later averaging.
pout=zeros(NumberAverages*length(p1),1);
yout=zeros(NumberAverages*length(y),1);
for n=1:NumberAverages
    pout((n-1)*N+1:n*N)=p1+max(p1)*NSRForce*randn(N,1);
    pout((n-1)*N+1:n*N)=pout((n-1)*N+1:n*N).*aforcew(N,10);
    yout((n-1)*N+1:n*N)=y+NSRResponse*max(abs(y))*randn(length(y),1);
end
[Gxx,Gyx,Gyy,f]=time2xmtrx(pout,yout,fs,boxcar(N),0);
[Ha,C]=xmtrx2frf(Gxx,Gyx,Gyy);
figure
subplot(2,2,1)
semilogy(f,abs(Ha),LineType{1},'LineWidth',LineWidth)
axis([fmin fmax amin amax])
grid
ylabel('Dyn. flexibility [m/N]','FontName',FontName,'FontSize',FontSize)
title('FRF without exp. win.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,3)
plot(f,C,LineType{1},'LineWidth',LineWidth)
axis([fmin fmax 0 1.1])
grid
ylabel('Coherence','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
title('No exp. win.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
% Apply exponential window
% Also to the force!
for n=1:NumberAverages
    pout((n-1)*N+1:n*N)=pout((n-1)*N+1:n*N).*aexpw(N,ExpWinEnd);
    yout((n-1)*N+1:n*N)=yout((n-1)*N+1:n*N).*aexpw(N,ExpWinEnd);
end
[Gxx,Gyx,Gyy,f]=time2xmtrx(pout,yout,fs,boxcar(N),0);
[Ha,C]=xmtrx2frf(Gxx,Gyx,Gyy);
subplot(2,2,2)
semilogy(f,abs(Ha),LineType{1},'LineWidth',LineWidth)
grid
ylabel('Dyn. flexibility [m/N]','FontName',FontName,'FontSize',FontSize)
axis([fmin fmax amin amax])
title('FRF with exp. win.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)
subplot(2,2,4)
plot(f,C,LineType{1},'LineWidth',LineWidth)
grid
axis([fmin fmax 0 1.1])
ylabel('Coherence','FontName',FontName,'FontSize',FontSize)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
title('With exp. win.','FontName',FontName,'FontSize',FontSize)
set(gca,'FontName',FontName,'FontSize',FontSize)

