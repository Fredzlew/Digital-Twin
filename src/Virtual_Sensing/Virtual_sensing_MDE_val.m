%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('data_sens'))

% Load data
% Low damping
% data = readmatrix('data_5_6_1.txt')'; % Loading displacement data
% High damping
% data = readmatrix('data_5_2_1.txt')'; % Loading displacement data
% filename = load('Eigenvalue_modeshape_residual_stiffmass.mat');
% U = filename.U;
% fss = data(2:6,:)/1000; % Converting mm to m
% xbase = data(7,:)/1000;
% xm = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]-xbase; % Swap columns due to sensor
% Simulated data
file = load('..\data\1_data_sim_newmark_jan.mat');
xm = file.dis_new;
U = file.U;

% Manually create timesteps (stepsize dt=0.001)
data = linspace(0,size(xm,2)*0.001,size(xm,2));

% Modeshapes from mottershead for High damping data
% U = [0.2338    0.6457    1.0000    1.0000   -0.6471;...
%     0.5108    1.0000    0.4949   -0.6345    1.0000;...
%     0.7398    0.6738   -0.8045   -0.4821   -0.9794;...
%     0.9066   -0.0908   -0.7113    0.9791    0.6360;...
%     1.0000   -0.7609    0.6006   -0.3793   -0.1715];



%% Virtual sensing part
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = 1;

% Index of measured locations (1 = bottom, 5 = top)
im = [5];

% Calculate displacement at predicted locations
[xp,qt] = VirtualSensVal(xm,U,num_ms,im);

% Number of time steps to plot
nt = 10000;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(xm(1,:)))
    plot(data(1,1:nt),xp(vs,1:nt),'r')
    plot(data(1,1:nt),xm(vs,1:nt),'b')
    legend('Predicted displacements','Actual displacements','FontSize',14)
else
    plot(data(1,1:nt),xm(vs,1:nt),'b')
    plot(data(1,1:nt),xp(vs,1:nt),'r')
    legend('Actual displacements','Predicted displacements','FontSize',14)
end
hold off
title(['Displacements at virtual sensor:',num2str(vs)],'FontSize',20)
xlabel('Time [s]','FontSize',14)
ylabel('Displacement [m]','FontSize',14)

%%% Calculate quality indicators %%%
% TRAC
TRAC = (xm(vs,:)*xp(vs,:)')^2/((xm(vs,:)*xm(vs,:)')*(xp(vs,:)*xp(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(TRAC)])

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(xm(vs,:)-xp(vs,:)))/size(xm,2)/std(xp(vs,:));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])

% Plot the modal coordinates in the frequency domain
Fs = 1000;           % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(data,2);    % Length of signal
t = (0:L-1)*T;       % Time vector
f = Fs*(0:(L/2))/L;

% Plots
Y = fft(qt(1,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")