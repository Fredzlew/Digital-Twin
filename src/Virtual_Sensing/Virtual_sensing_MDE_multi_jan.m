%% Script that performs virtual sensing %%
clear;clc;close all

% % Adding path to data
% addpath(genpath('data_sens'))
% 
% % Load data
% file = load('.\data\simulated_data.mat');
% xm = file.xmnew;
% filename = load('.\data\modelprop_jan.mat');
% U = filename.U;
% 
% % Manually create timesteps (stepsize dt=0.001)
% data = linspace(0,size(xm,2)*0.001,size(xm,2));

%.......... STRUCTURE .....................................................
%.......... STRUCTURE .....................................................
rng(1);
% dofs (floors)
n = 5;

% discrete parameters
m  = 1;                                 % (lumped) mass
k  = m*pi^2/sin(pi/2/(2*n+1))^2;        % stiffness

% stiffness matrix
KK      = 2*diag(ones(n,1)) - diag(ones(n-1,1),-1)...
         - diag(ones(n-1,1),1);
K = k*KK;
K(n,n)    = k;                          % free end at top floor

% mass matrix
M = m*eye(n);

% damping matrix
C = 0.00*K;

% solving eigenvalue problem
[V,D] = eig(K,M);
wj    = sqrt(diag(D))
freq_f = wj/(2*pi)

% normalized 
for jj=1:n
   [~,imax] = max(abs(V(:,jj)));
   V(:,jj) = V(:,jj)/V(imax,jj);
   %V(:,jj) = V(:,jj)/V(end,jj);
end

% load
dt = 0.001;
N = 50000;

% load
%ff = ones(n,1)*randn(1,N+1);
ff = [ 0 ; 0 ; 0 ; 0 ; 1 ]*randn(1,N+1);

x0 = zeros(n,1);
v0 = x0;

[xm,vm,am,t] = newmark(K,C,M,x0,v0,dt,N,ff,0.25,0.5);

figure(1)
plot(t,xm(1,:),'b.-')

x1 = xm(1,:);
fs = 1/dt;
f_cut = 4;

x1_0_cut = bandpass(x1,[0.5 f_cut],fs);

x1_cut_100 = bandpass(x1,[f_cut 100],fs);

x1_fil = x1_0_cut + x1_cut_100;

hold on
plot(t,x1_fil,'r.-')

hold off
% TRAC
TRAC = (x1*x1_fil')^2/((x1*x1')*(x1_fil*x1_fil'));
disp(['TRAC value for sensor ',num2str(1),' with ',num2str(0),' modes:',num2str(TRAC)])

xm_0_cut = bandpass(xm,[0.5 f_cut],fs);

xm_cut_100 = bandpass(xm,[f_cut 100],fs);
%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1,2,3];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Calculate displacement at predicted locations
[xp1,qt1] = VirtualSensVal(xm_0_cut,V,num_ms,im);

% Virtual sensing part 2
% Calculate displacement at predicted locations
[xp2,qt2] = VirtualSensVal(xm_cut_100,V,num_ms,im);

xp = xp1 + xp2;
qt = qt1 + qt2;

%% Plotting
% Number of time steps to plot
nt = 10000;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(x1))
    plot(t(1:nt),xp(vs,1:nt),'r')
    plot(t(1:nt),x1(1:nt),'b')
    legend('Predicted displacements','Actual displacements','FontSize',14)
else
    plot(t(1:nt),x1(1:nt),'b')
    plot(t(1:nt),xp(vs,1:nt),'r')
    legend('Actual displacements','Predicted displacements','FontSize',14)
end
hold off
title(['Displacements at virtual sensor:',num2str(vs)],'FontSize',20)
xlabel('Time [s]','FontSize',14)
ylabel('Displacement [m]','FontSize',14)

%%% Calculate quality indicators %%%
% TRAC
TRAC = (x1*xp(vs,:)')^2/((x1*x1')*(xp(vs,:)*xp(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(TRAC)])

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(x1-xp(vs,:)))/size(x1,2)/std(xp(vs,:));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])

% Plot the modal coordinates in the frequency domain
Fs = 1000;           % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(t,2);    % Length of signal
f = Fs*(0:(L/2))/L;

% Plots
Y = fft(qt(3,:)); % mode vi plotter
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")