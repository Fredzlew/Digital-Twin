%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('..\..\data'),genpath('..\..\npy-matlab-master'),genpath('.\Filtered_data'),genpath('..\function'),genpath('..\..\Model_updating\Sensitivity_method\data_updated_par_sens'))

promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);
if q == 1
    % Load data
    xm_mode1 = readNPY('.\Filtered_data\data_filt_mode1_high.npy');
    xm_mode2 = readNPY('.\Filtered_data\data_filt_mode2_high.npy');
    xm_mode3 = readNPY('.\Filtered_data\data_filt_mode3_high.npy');
    xm_mode4 = readNPY('.\Filtered_data\data_filt_mode4_high.npy');
    xm_mode5 = readNPY('.\Filtered_data\data_filt_mode5_high.npy');
    xm_filt = readNPY('.\Filtered_data\data_filt_all_high.npy');
    xm_data =  xm_mode1+xm_mode2+xm_mode3+xm_mode4+xm_mode5;
    % loading mode shapes
    promptt = "1 = residual Eig + MS, 2 = residual Eig, 3 = resiudal MS, 4 = Eig, 5 = MS, 6 = MS (MAC), 7 = Eig + MS, 8 = Eig + MS (EIL): ";
    zz = input(promptt);
    if zz == 1
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_Mode_shape_residual_high.mat');
    elseif zz == 2
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_residual_high.mat');
    elseif zz == 3
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Mode_shape_residual_high.mat');
    elseif zz == 4
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreq_high.mat');
    elseif zz == 5
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSImode_high.mat');
    elseif zz == 6
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSImode_mac_high.mat');
    elseif zz == 7
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreqmode_high.mat');
    elseif zz == 8
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreqmodeEIL_high.mat');
    end
    U = filename.U;
elseif q == 2
    % Load data
    xm_mode1 = readNPY('.\Filtered_data\data_filt_mode1_no_damp.npy');
    xm_mode2 = readNPY('.\Filtered_data\data_filt_mode2_no_damp.npy');
    xm_mode3 = readNPY('.\Filtered_data\data_filt_mode3_no_damp.npy');
    xm_mode4 = readNPY('.\Filtered_data\data_filt_mode4_no_damp.npy');
    xm_mode5 = readNPY('.\Filtered_data\data_filt_mode5_no_damp.npy');
    xm_filt = readNPY('.\Filtered_data\data_filt_all_no_damp.npy');
    xm_data =  xm_mode1+xm_mode2+xm_mode3+xm_mode4+xm_mode5;
    % loading mode shapes
    promptt = "1 = residual Eig + MS, 2 = residual Eig, 3 = resiudal MS, 4 = Eig, 5 = MS, 6 = MS (MAC), 7 = Eig + MS, 8 = Eig + MS (EIL): ";
    zz = input(promptt);
    if zz == 1
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat');
    elseif zz == 2
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_residual_no_damp.mat');
    elseif zz == 3
        filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Mode_shape_residual_no_damp.mat');
    elseif zz == 4
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreq_no_damp.mat');
    elseif zz == 5
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSImode_no_damp.mat');
    elseif zz == 6
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSImode_mac_no_damp.mat');
    elseif zz == 7
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreqmode_no_damp.mat');
    elseif zz == 8
        filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreqmodeEIL_no_damp.mat');
    end
   
    U = filename.U;
end

% Make a time vector
Fs = 50;             % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(xm_filt,2);  % Length of signal
t = (0:L-1)*T;       % Time vector

% filename = load('..\..\data\modelprop.mat');
% U = filename.U;

% Verifying original and filtered data

% Number of time steps to plot
nt = 500;
% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

hold on 
%plot(t(1:nt),xm_data(vs,1:nt),'b.-')
%plot(t(1:nt),xm_filt(vs,1:nt),'r.-')
hold off
% TRAC
TRAC = (xm_data(vs,:)*xm_filt(vs,:)')^2/((xm_data(vs,:)*xm_data(vs,:)')*(xm_filt(vs,:)*xm_filt(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(0),' modes:',num2str(TRAC)])
%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1];

% Index of measured locations (1 = bottom, 5 = top)
im = [3,5];

% Calculate displacement at predicted locations
[xp1,qt1] = VirtualSensVal(xm_mode1,U,num_ms,im);

% Virtual sensing part 2
% Number of modeshapes included in approximation (max length(im))
num_ms = [2];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,5];
% Calculate displacement at predicted locations
[xp2,qt2] = VirtualSensVal(xm_mode2,U,num_ms,im);

% Virtual sensing part 3
% Number of modeshapes included in approximation (max length(im))
num_ms = [3];

% Index of measured locations (1 = bottom, 5 = top)
im = [3,5];
% Calculate displacement at predicted locations
[xp3,qt3] = VirtualSensVal(xm_mode3,U,num_ms,im);


% Virtual sensing part 4
% Number of modeshapes included in approximation (max length(im))
num_ms = [4];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,5];
% Calculate displacement at predicted locations
[xp4,qt4] = VirtualSensVal(xm_mode4,U,num_ms,im);

% Virtual sensing part 5
% Number of modeshapes included in approximation (max length(im))
num_ms = [5];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3];
% Calculate displacement at predicted locations
[xp5,qt5] = VirtualSensVal(xm_mode5,U,num_ms,im);

xp = xp1 + xp2 + xp3 + xp4 + xp5;
%qt = qt1 + qt2;

%% Plotting
% Number of time steps to plot
nt = 500;

% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

% Plot actual displacements vs predicted displacements
figure
hold on
if sum(abs(xp(1,:)))>sum(abs(xm_filt(1,:)))
    plot(t(1:nt),xp(vs,1:nt),'r')
    plot(t(1:nt),xm_filt(vs,1:nt),'b')
    legend('Predicted displacements','Actual displacements','FontSize',14)
else
    plot(t(1:nt),xm_filt(vs,1:nt),'b')
    plot(t(1:nt),xp(vs,1:nt),'r')
    legend('Actual displacements','Predicted displacements','FontSize',14)
end
hold off
title(['Displacements at virtual sensor:',num2str(vs)],'FontSize',20)
xlabel('Time [s]','FontSize',14)
ylabel('Displacement [m]','FontSize',14)

%%% Calculate quality indicators %%%
% TRAC
TRAC = (xm_filt(vs,:)*xp(vs,:)')^2/((xm_filt(vs,:)*xm_filt(vs,:)')*(xp(vs,:)*xp(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(TRAC)])

% FRAC
FRAC = (fft(xm_filt(vs,:))*ctranspose(fft(xp(vs,:))))^2/((fft(xm_filt(vs,:))*ctranspose(fft(xm_filt(vs,:))))*(fft(xp(vs,:))*ctranspose(fft(xp(vs,:)))));
disp(['FRAC value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(FRAC)])

% MAE (normalized with respect to the standard deviation)
MAE = sum(abs(xm_filt(vs,:)-xp(vs,:)))/size(xm_filt,2)/std(xp(vs,:));
disp(['MAE value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(MAE)])

% ME (mean error)
ME = sum(xm_filt(vs,:)-xp(vs,:))/size(xm_filt,2);
disp(['ME value for sensor ',num2str(vs),' with ',num2str(num_ms),' modes:',num2str(ME)]) 

% Plot the modal coordinates in the frequency domain
f = Fs*(0:(L/2))/L;

% Plots
Y = fft(xm_filt(1,:)); % mode vi plotter
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")

% Plots
Y = fft(xp(1,:)); % mode vi plotter
P2 = abs(Y/L);
P11 = P2(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

figure
plot(f,P11) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 20])
ylabel("|P1(f)|")


% Plots meaasered - ori data
xx = xp(1,:)-xm_filt(1,:);
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
%% Download the file to plot in latex
%{
if q == 1
    T1 = array2table([num2cell(xp(1,:)'),num2cell(xm_filt(1,:)')]);
    T1.Properties.VariableNames(1:2) = {'x_pred','x_meas'};
    writetable(T1,'.\Filtered_data\data_predicted_high.txt','Delimiter',' ')
elseif q == 2
    T1 = array2table([num2cell(xp(1,:)'),num2cell(xm_filt(1,:)')]);
    T1.Properties.VariableNames(1:2) = {'x_pred','x_meas'};
    writetable(T1,'.\Filtered_data\data_predicted_nodamp.txt','Delimiter',' ')
end
%}
