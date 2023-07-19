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
    xm_data =  readNPY('.\Filtered_data\data_filt_high.npy')';
    % loading mode shapes
    filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_Mode_shape_residual_high.mat');
    U = filename.U;
elseif q == 2
    % Load data
    xm_mode1 = readNPY('.\Filtered_data\data_filt_mode1_no_damp.npy');
    xm_mode2 = readNPY('.\Filtered_data\data_filt_mode2_no_damp.npy');
    xm_mode3 = readNPY('.\Filtered_data\data_filt_mode3_no_damp.npy');
    xm_mode4 = readNPY('.\Filtered_data\data_filt_mode4_no_damp.npy');
    xm_mode5 = readNPY('.\Filtered_data\data_filt_mode5_no_damp.npy');
    xm_filt = readNPY('.\Filtered_data\data_filt_all_no_damp.npy');
    xm_data =  readNPY('.\Filtered_data\data_filt_no_damp.npy')';
    % loading mode shapes
    filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat');
    %filename = load('..\..\Model_updating\Costfunctions\data_updated_par\SSIfreqmode_no_damp.mat');
    %filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Mode_shape_residual_no_damp.mat');
    %filename = load('..\..\Model_updating\Sensitivity_method\data_updated_par_sens\Eigenvalue_residual_no_damp.mat');
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
plot(t(1:nt),xm_data(vs,1:nt),'b.-')
plot(t(1:nt),xm_filt(vs,1:nt),'r.-')
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

%%
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
%% save the predicted displacement as txt
%{
if q == 1
    T1 = array2table(num2cell(xp(1,:)'));
    T1.Properties.VariableNames(1) = {'x_pred'};
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src\Virtual_sensing\Experimental\Filtered_data\data_predicted_high.txt','Delimiter',' ')
elseif q == 2
    T1 = array2table(num2cell(xp(1,:)'));
    T1.Properties.VariableNames(1) = {'x_pred'};
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src\Virtual_sensing\Experimental\Filtered_data\data_predicted_nodamp.txt','Delimiter',' ')
end
%}

%% download data
if q == 1
    T_pred = array2table([num2cell(t(1:500)'),num2cell(xp(1,1:500)'),num2cell(xm_filt(1,1:500)'),num2cell(xm_data(1,1:500)')]);
    T_fft = array2table([num2cell(decimate(f(1:604801)',60)),num2cell(decimate(P1(1:604801)',60)),num2cell(decimate(P11(1:604801)',60)),num2cell(decimate(P1_pred_meas(1:604801)',60))]);
    
    T_pred.Properties.VariableNames(1:4) = {'data','pred','meas','ori'};
    T_fft.Properties.VariableNames(1:4) = {'f','P1_meas','P1_pred','P1_pred_meas'};
    
    
    writetable(T_pred,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_pred_multi5_highdamp.csv','Delimiter',';')
    writetable(T_fft,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_fft_multi5_highdamp.csv','Delimiter',';')
elseif q == 2
    T_pred = array2table([num2cell(t(1:500)'),num2cell(xp(1,1:500)'),num2cell(xm_filt(1,1:500)'),num2cell(xm_data(1,1:500)')]);
    T_fft = array2table([num2cell(decimate(f(1:604801)',60)),num2cell(decimate(P1(1:604801)',60)),num2cell(decimate(P11(1:604801)',60)),num2cell(decimate(P1_pred_meas(1:604801)',60))]);
    
    T_pred.Properties.VariableNames(1:4) = {'data','pred','meas','ori'};
    T_fft.Properties.VariableNames(1:4) = {'f','P1_meas','P1_pred','P1_pred_meas'};
    
    
    writetable(T_pred,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_pred_multi5_nodamp.csv','Delimiter',';')
    writetable(T_fft,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_fft_multi5_nodamp.csv','Delimiter',';')
end