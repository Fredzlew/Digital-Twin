%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('data_sens'))

% Load data
% No damping
data = readmatrix('..\data\Anela\data.txt'); % Loading displacement data
data = downsample(data,2)';
xm_0_cut = readNPY('..\python\data\Modal_parameters_anela\data_filt_0_cut.npy')'/1000;
xm_cut_end = readNPY('..\python\data\Modal_parameters_anela\data_filt_cut_end.npy')'/1000;

xm_ori = readNPY('..\python\data\Modal_parameters_anela\data_nodamp_filtdata.npy')'/1000;


% data = readmatrix('data_5_2_1.txt')'; % Loading displacement data
filename = load('data_sens\Eigenvalue_modeshape_residual_stiff_nodamp.mat');
U = filename.U;

% Verifying original and filtered data
xm_filt =  xm_0_cut + xm_cut_end;

% Number of time steps to plot
nt = 500;
% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

hold on 
plot(data(1,1:nt),xm_ori(vs,1:nt),'b.-')
plot(data(1,1:nt),xm_filt(vs,1:nt),'r.-')
hold off
% TRAC
TRAC = (xm_ori(vs,:)*xm_filt(vs,:)')^2/((xm_ori(vs,:)*xm_ori(vs,:)')*(xm_filt(vs,:)*xm_filt(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(0),' modes:',num2str(TRAC)])
%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1,2,3];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Calculate displacement at predicted locations
[xp1,qt1] = VirtualSensVal(xm_0_cut,U,num_ms,im);

% Virtual sensing part 2
% Calculate displacement at predicted locations
[xp2,qt2] = VirtualSensVal(xm_cut_end,U,num_ms,im);

xp = xp1 + xp2;
qt = qt1 + qt2;

%% Plotting
% Number of time steps to plot
nt = 500;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(xm_ori(1,:)))
    plot(data(1,1:nt),xp(vs,1:nt),'r')
    plot(data(1,1:nt),xm_ori(vs,1:nt),'b')
    legend('Predicted displacements','Actual displacements','FontSize',14)
else
    plot(data(1,1:nt),xm_ori(vs,1:nt),'b')
    plot(data(1,1:nt),xp(vs,1:nt),'r')
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

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(xm_ori(vs,:)-xp(vs,:)))/size(xm_ori,2)/std(xp(vs,:));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])

% Plot the modal coordinates in the frequency domain
Fs = 50;           % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(data,2);    % Length of signal
t = (0:L-1)*T;       % Time vector
f = Fs*(0:(L/2))/L;

% Plots
Y = fft(qt(1,:)); % mode vi plotter
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")