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
n=size(f,1);
dt=1/fs; %sampling rate
% Solve eigenvalue problem to find numerical modal parameters
[Us, Values]=eig(K,M);
Freq=sqrt(diag(Values))/(2*pi); % undamped natural frequency
steps=size(f,2);
omegas = filename.fn;

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
for i=100:10:1000 % Doesn't find 5 eigenvalues above 2600
    nm = 5; %Number of modes
    Y=f; %Displacements
    nrows=i;%112;%50*(2*nm/5)+1;     %more than 20 * number of modes
    ncols=7400;%7044;%4/5*size(f,2)-nrows-3;    %more than 2/3 of No. of data
    inputs=1;     
    cut=2*nm;        %Identify 5 modes
    shift=10;      %Adjust EMAC sensitivity
    EMAC_option=1; %EMAC is calculated only from observability matrix

    [Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA

    OMAfreq(:,i/100) = Result.Parameters.NaFreq;
    Acc(i/100) = min(sum(OMAfreq),sum(omegas))/max(sum(OMAfreq),sum(omegas));
end
figure
plot(linspace(1,100,100),Acc)
title('Convergence analysis of ERA')
xlabel('Iteration')
ylabel('Accuracy')