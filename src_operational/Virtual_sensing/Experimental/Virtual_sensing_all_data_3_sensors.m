%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('..\..\data'),genpath('.\Filtered_data_3_sensors'),genpath('..\function'),genpath('..\..\Model_updating\Sensitivity_method\data_updated_par_sens_3_sensors'))

promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);
if q == 1
    % Load data
    xm_0_cut = readNPY('.\Filtered_data_3_sensors\data_filt_all_3_sensors_high.npy');
    %xm_ori = readNPY('.\Filtered_data_3_sensors\data_filt_all_3_sensors_high.npy');
    xm_ori = readNPY('.\Filtered_data_3_sensors\data_filt_all_sensors_high.npy');
    % loading mode shapes
    filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens_3_sensors\Eigenvalue_Mode_shape_residual_3_sensors_high.mat');
    U = filename.U;
elseif q == 2
    % Load data
    xm_0_cut = readNPY('.\Filtered_data_3_sensors\data_filt_all_3_sensors_no_damp.npy');
    %xm_ori = readNPY('.\Filtered_data_3_sensors\data_filt_all_3_sensors_no_damp.npy');
    xm_ori = readNPY('.\Filtered_data_3_sensors\data_filt_all_sensors_no_damp.npy');
    % loading mode shapes
    filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens_3_sensors\Eigenvalue_Mode_shape_residual_3_sensors_no_damp.mat');
    U = filename.U;
end

% Make a time vector
Fs = 50;             % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(xm_ori,2);  % Length of signal
t = (0:L-1)*T;       % Time vector

% filename = load('..\..\data\modelprop.mat');
% U = filename.U;

%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1,2];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Which sensor (1 = sensor 2, 2 = sensor 3, 3 = sensor 5)
ls = [1,2,3];

% Calculate displacement at predicted locations
[xp,qt] = VirtualSensVal_3_sensors(xm_0_cut,U,num_ms,im,ls);


%% Plotting
% Number of time steps to plot
nt = 500;

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

% ME (mean error)
ME = sum(xm_ori(vs,:)-xp(vs,:))/size(xm_ori,2);
disp(['ME value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(ME)]) 

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
xlim([0 14])
ylabel("|P1(f)|")

% Plots
Y = fft(xp(1,:)); % mode vi plotter
P22 = abs(Y/L);
P11 = P22(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

figure
plot(f,P11) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 14])
ylabel("|P1(f)|")

% Plots pred - ori
xx = xp(1,:)-xm_ori(1,:);
Y = fft(xx); % mode vi plotter
P2 = abs(Y/L);
P1_pred_meas = P2(1:L/2+1);
P1_pred_meas(2:end-1) = 2*P1_pred_meas(2:end-1);

figure
plot(f,P1_pred_meas) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")