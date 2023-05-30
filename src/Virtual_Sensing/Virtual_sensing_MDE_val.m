%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('data_sens'))

% Load data
% High damping
data = readmatrix('data_5_2_1.txt')'; % Loading displacement data
% Low damping
% data = readmatrix('data_5_6_1.txt')'; % Loading displacement data

fss = data(2:6,:)/1000; % Converting mm to m
xm = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
filename = load('Eigenvalue_modeshape_residual_stiffmass.mat');
U = filename.U;
% U = [0.3164, 0.7748, 1.0000, 1.0000, -0.4440;...
%     0.6301, 1.0000, 0.2081, -0.9371,0.8281;...
%     0.7783, 0.4530, -0.7971, -0.0186, -0.9820;...
%     1.0000, -0.3185, -0.3893, 0.8750, 1.0000;...
%     0.9923, -0.7864, 0.6152, -0.5075, -0.3861];



%% Virtual sensing part
%close all;
% Number of modeshapes included in approximation
num_ms = 3;

% Index of measured locations (1 = bottom, 5 = top)
im = [3,4,5];

% Calculate displacement at predicted locations
[xp] = VirtualSensVal(xm,U,num_ms,im);

% Number of time steps to plot
nt = 4320000;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(xm(1,:)))
    plot(data(1,1:nt),xp(vs,1:nt),'r')
    plot(data(1,1:nt),xm(vs,1:nt),'b')
    legend('Predicted displacements','Actual displacements')
else
    plot(data(1,1:nt),xm(vs,1:nt),'b')
    plot(data(1,1:nt),xp(vs,1:nt),'r')
    legend('Actual displacements','Predicted displacements')
end
hold off
title('Displacements at virtual sensor')
xlabel('Time [s]')
ylabel('Displacement [m]')

%%% Calculate quality indicators %%%
% TRAC
TRAC = (xm(vs,1:nt)*xp(vs,1:nt)')^2/(xm(vs,1:nt)*xm(vs,1:nt)'*xp(vs,1:nt)*xp(vs,1:nt)');
disp(['TRAC value: ',num2str(TRAC)])

% MAE (in meters)
MAE = sum(abs(xm(vs,1:nt)-xp(vs,1:nt)))/nt;
disp(['MAE value: ',num2str(MAE)])
