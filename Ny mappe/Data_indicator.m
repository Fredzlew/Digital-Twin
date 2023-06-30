%% Script that investigate each sensor %%
clear;clc;close all
data = readmatrix('data\Anela\data.txt')'; % Loading displacement data
fss = data(2:6,:); % [mm]
xm = flip(fss,1); % Swap rows due to sensor


xm_std = zeros(5,1);
xm_mean = zeros(5,1);

for i = 1:5
xm_std(i) = std(xm(i,:));
xm_mean(i) = mean(abs(xm(i,:)));
end

% Using filter data
data = downsample(data',2)';
xm = readNPY('.\python\data\Modal_parameters_anela\data_nodamp_filtdata.npy')'/1000;

%% FFT of the data

% modal coordinates
phi_m = readNPY('.\python\data\Modal_parameters_anela\SSIphi_nodamp.npy');
% Calculate the psuedo-inverse of the measured mode shapes
phi_minv = (phi_m'*phi_m)^-1*phi_m';
% phi_minv = phi_m'*inv(phi_m*phi_m')
% phi_minv = pinv(phi_m)

% Calculate modal coordinates
qt = phi_minv*xm;


% Plot the modal coordinates in the frequency domain
Fs = 50;            % Sampling frequency                    
T = 1/Fs;            % Sampling period       
L = size(data,2);    % Length of signal
t = data(1,:);       % Time vector
f = Fs*(0:(L/2))/L; % frequency domain

% Plots
Y = fft(qt(1,:)); % mode vi plotter
P2 = abs(Y/L); % two-sided spectrum
P1 = P2(1:L/2+1); % single-sided amplitude spectrum P1
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
xlim([0 15])
ylabel("|P1(f)|")