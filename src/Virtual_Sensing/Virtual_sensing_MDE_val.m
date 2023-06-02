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
% U = [0.2338    0.6457    1.0000    1.0000   -0.6471;...
%     0.5108    1.0000    0.4949   -0.6345    1.0000;...
%     0.7398    0.6738   -0.8045   -0.4821   -0.9794;...
%     0.9066   -0.0908   -0.7113    0.9791    0.6360;...
%     1.0000   -0.7609    0.6006   -0.3793   -0.1715];



%% Virtual sensing part
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = 2;

% Index of measured locations (1 = bottom, 5 = top)
im = [4,5];

% Calculate displacement at predicted locations
[xp] = VirtualSensVal(xm,U,num_ms,im);

% Number of time steps to plot
nt = 100;

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
title(['Displacements at virtual sensor:',num2str(vs)])
xlabel('Time [s]')
ylabel('Displacement [m]')

%%% Calculate quality indicators %%%
% TRAC
TRAC = (xm(vs,1:nt)*xp(vs,1:nt)')^2/((xm(vs,1:nt)*xm(vs,1:nt)')*(xp(vs,1:nt)*xp(vs,1:nt)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(TRAC)])

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(xm(vs,1:nt)-xp(vs,1:nt)))/nt/std(xp(vs,1:nt));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])
