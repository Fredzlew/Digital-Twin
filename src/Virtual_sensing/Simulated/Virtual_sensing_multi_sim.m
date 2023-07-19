%% Script that performs virtual sensing %%
clear;clc;close all

% Adding path to data
addpath(genpath('..\..\data'),genpath('.\Filtered_data_sim'),genpath('..\function'))

% Load data
xm_0_cut = readNPY('.\Filtered_data_sim\data_filt_0_cut_sim.npy');
xm_cut_end = readNPY('.\Filtered_data_sim\data_filt_cut_end_sim.npy');

xm_ori = readNPY('.\Filtered_data_sim\data_filtdata_all_sim.npy');

% Make a time vector
Fs = 500;            % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(xm_ori,2);  % Length of signal
t = (0:L-1)*T;       % Time vector


% data = readmatrix('data_5_2_1.txt')'; % Loading displacement data
filename = load('..\..\data\modelprop.mat');
U = filename.U;

% Verifying original and filtered data
xm_filt =  xm_0_cut + xm_cut_end;

% Number of time steps to plot
nt = 1000;
% Show displacements for # virtual sensor (1 = bottom)
vs = 1;

hold on 
plot(t(1:nt),xm_ori(vs,1:nt),'b.-')
plot(t(1:nt),xm_filt(vs,1:nt),'r.-')
hold off
% TRAC
TRAC = (xm_ori(vs,:)*xm_filt(vs,:)')^2/((xm_ori(vs,:)*xm_ori(vs,:)')*(xm_filt(vs,:)*xm_filt(vs,:)'));
disp(['TRAC value for sensor ',num2str(vs),' with ',num2str(0),' modes:',num2str(TRAC)])
%% Virtual sensing part 1
%close all;
% Number of modeshapes included in approximation (max length(im))
num_ms = [1,2];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Calculate displacement at predicted locations
[xp1,qt1] = VirtualSensVal(xm_0_cut,U,num_ms,im);

% Virtual sensing part 2
% Number of modeshapes included in approximation (max length(im))
num_ms = [3,4];

% Index of measured locations (1 = bottom, 5 = top)
im = [2,3,5];

% Calculate displacement at predicted locations
[xp2,qt2] = VirtualSensVal(xm_cut_end,U,num_ms,im);

xp = xp1 + xp2;
%qt = qt1 + qt2;

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
%% download data
T_pred = array2table([num2cell(t(1:5000)'),num2cell(xp(1,1:5000)'),num2cell(xm_ori(1,1:5000)')]);
T_fft = array2table([num2cell(decimate(f(1:180001)',17)),num2cell(decimate(P1(1:180001)',17)),num2cell(decimate(P11(1:180001)',17)),num2cell(decimate(P1_pred_meas(1:180001)',17))]);

T_pred.Properties.VariableNames(1:3) = {'data','pred','meas'};
T_fft.Properties.VariableNames(1:4) = {'f','P1_meas','P1_pred','P1_pred_meas'};


writetable(T_pred,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_pred_multi_sim.csv','Delimiter',';')
writetable(T_fft,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_fft_multi_sim.csv','Delimiter',';')