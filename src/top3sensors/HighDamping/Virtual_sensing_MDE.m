%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('data_sens'))

% Load data
data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
fss = data(2:6,:)/1000; % Converting mm to m
f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
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
num_ms = 2;

% Index of measured locations (1 = bottom, 5 = top)
im = [3,4,5];

% Calculate displacement at predicted locations
[xp] = VirtualSensVal(f,U,num_ms,im);

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(2,:)))>sum(abs(f(2,:)))
    plot(data(1,1:60000),xp(2,1:60000),'r')
    plot(data(1,1:60000),f(2,1:60000),'b')
    legend('Predicted displacements','Actual displacements')
else
    plot(data(1,1:60000),f(2,1:60000),'b')
    plot(data(1,1:60000),xp(2,1:60000),'r')
    legend('Actual displacements','Predicted displacements')
end
hold off
title('Displacements at first virtual sensor')
xlabel('Time [s]')
ylabel('Displacement [m]')

% Plot the relative error between actual displacements and predicted displacements
% err = zeros(size(xp,1),size(xp,2));
% for j = 1:size(xp,1)
%     for i = 1:size(xp,2)
%         err(j,i) = abs(abs(xp(j,i))-abs(f(j,i)))/abs(f(j,i))*100;
%     end
% end
% 
% figure (2)
% plot(err(1,:))
% title('Relative error')
% xlabel('Time step [-]')
% ylabel('Error [%]')









