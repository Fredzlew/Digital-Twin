clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Choose of data
prompt = "Use ERA for measured or simulated data (1=measured, 2=simulated)? ";
ERAdata = input(prompt);
if ERAdata == 1
    % Measurements
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    fss = data(2:6,1:10000)/1000; % Converting mm to m
    f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
elseif ERAdata == 2
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
end

filename = load('modelprop.mat'); % Loads mass and stiffness matrices
omega_min = 1.70; % Minimum natural frequency
zeta_min = 0.015; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=0;%alpha*M+beta*K; % Damping matrix using Rayleigh damping
fs=100; % Sampling frequency (1/dt)
omegas = filename.fn;
%Apply modal superposition to get response
%--------------------------------------------------------------------------

n=size(f,1); % Number of floors/sensors
dt=1/fs; %sampling rate
[Us, Values]=eig(K,M); % Solving eigenvalueproblem (undamped)
Freq=sqrt(diag(Values))/(2*pi); % Undamped natural frequency
steps=size(f,2); % Number of samples

% normalizing mode shapes
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        Vectors(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization



%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
y=0;
for i=100:10:4000
    nm = 5; %number of modes
    output=f; % Displacements
    ncols=7400;%4/5*length(f); % More than 2/3*number of samples
    nrows=i;%50*nm; % More than 20*number of sensors
    cut=2*nm;  % cut=4 -> 2 modes, cut=10 -> 5 modes
    [Result]=SSID(output,fs,ncols,nrows,cut);    %SSI
    y=y+1;
    iv(y) = i;
    OMAfreq(:,y) = Result.Parameters.NaFreq;
    Diff(:,y) = [min(OMAfreq(1,y),omegas(1))/max(OMAfreq(1,y),omegas(1));min(OMAfreq(2,y),omegas(2))/max(OMAfreq(2,y),omegas(2));min(OMAfreq(3,y),omegas(3))/max(OMAfreq(3,y),omegas(3));min(OMAfreq(4,y),omegas(4))/max(OMAfreq(4,y),omegas(4));min(OMAfreq(5,y),omegas(5))/max(OMAfreq(5,y),omegas(5))]*100;
    Acc(y) = sum(Diff(:,y));
end
figure
plot(iv,Diff)
title('Convergence analysis of SSI')
xlabel('Number of rows in hankel matrix')
ylabel('Accuracy')
legend('f1','f2','f3','f4','f5')