%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('..\..\data'),genpath('.\Filtered_data_sim'),genpath('..\function'))

% Load data
xm_20 = readNPY('.\Filtered_data_sim\data_filtdata_all_sim.npy');
xm_ori = readNPY('.\Filtered_data_sim\sim_data.npy')';

% Make a time vector
Fs = 500;            % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(xm_ori,2);  % Length of signal
t = (0:L-1)*T;       % Time vector


% data = readmatrix('data_5_2_1.txt')'; % Loading displacement data
filename = load('..\..\data\modelprop.mat');
U = filename.U;

%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1,2];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Calculate displacement at predicted locations
[xp,qt] = VirtualSensVal(xm_20,U,num_ms,im);


%% Plotting
% Number of time steps to plot
nt = 10000;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(xm_ori(1,:)))
    plot(t(1:nt),xp(vs,1:nt),'r')
    plot(t(1:nt),xm_ori(vs,1:nt),'b')
    legend('Predicted displacements','Actual displacements','FontSize',14)
else
    plot(t(1:nt),xm_ori(vs,1:nt),'b')
    plot(t(1:nt),xp(vs,1:nt),'r')
    legend('Actual displacements','Predicted displacements','FontSize',14)
end
hold off
title(['Displacements at virtual sensor:',num2str(vs)],'FontSize',20)
xlabel('Time [s]','FontSize',14)
ylabel('Displacement [m]','FontSize',14)

%%% Calculate quality indicators %%%
% TRAC
TRAC = (xm_ori(vs,:)*xp(vs,:)')^2/((xm_ori(vs,:)*xm_ori(vs,:)')*(xp(vs,:)*xp(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(TRAC)])

% FRAC
FRAC = (fft(xm_ori(vs,:))*ctranspose(fft(xp(vs,:))))^2/((fft(xm_ori(vs,:))*ctranspose(fft(xm_ori(vs,:))))*(fft(xp(vs,:))*ctranspose(fft(xp(vs,:)))));
disp(['FRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(FRAC)])

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(xm_ori(vs,:)-xp(vs,:)))/size(xm_ori,2)/std(xp(vs,:));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])

% Plot the modal coordinates in the frequency domain
f = Fs*(0:(L/2))/L;

% Plots
Y = fft(xm_ori(1,:)); % mode vi plotter
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")
%% xm ori modal coordinates
% Calculate the psuedo-inverse of the measured mode shapes
U_inv = (U'*U)^-1*U';
% phi_minv = phi_m'*inv(phi_m*phi_m')
% phi_minv = pinv(phi_m)

% Calculate modal coordinates
qt = U_inv*xm_ori;
% Plots
Y = fft(qt(4,:)); % mode vi plotter
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")