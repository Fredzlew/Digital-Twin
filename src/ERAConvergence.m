clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Choose stiffness matrix
prompttt = "Ricky and Johan's or Jan's stiffness matrix: (1=Ricky&Johan, 2=Jan)? ";
prop = input(prompttt);
if prop == 1
    filename = load('modelprop.mat'); % Loads mass and stiffness matrices
elseif prop == 2
    filename = load('modelprop_jan.mat'); % Loads mass and stiffness matrices
end

% Choose of data
prompt = "Use ERA for measured or simulated data (1=measured, 2=simulated(IRF), 3=simulated(Newmark))? ";
ERAdata = input(prompt);
if ERAdata == 1
    % Measurements
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    fss = data(2:6,1:10000)/1000; % Converting mm to m
    f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
elseif ERAdata == 2 && prop == 1
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
elseif ERAdata == 2 && prop == 2
    % Simulated data
    data_sim = load('data_sim_jan.mat');
    f = data_sim.dis(:,1:10000);
elseif ERAdata == 3 && prop == 1
    % Simulated data
    data_sim = load('data_sim_newmark.mat');
    f = data_sim.dis_new(:,1:10000);
elseif ERAdata == 3 && prop == 2
    % Simulated data
    data_sim = load('data_sim_newmark_jan.mat');
    f = data_sim.dis_new(:,1:10000);
end

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
y=0;
for i=10:10:600 % Doesn't find 5 eigenvalues above 2600.
    nm = 5; %Number of modes
    Y=f; %Displacements
    nrows=i;%112;%50*(2*nm/5)+1;     %more than 20 * number of modes
    ncols=7400;%7044;%4/5*size(f,2)-nrows-3;    %more than 2/3 of No. of data
    inputs=1;     
    cut=2*nm;        %Identify 5 modes
    shift=10;      %Adjust EMAC sensitivity
    EMAC_option=1; %EMAC is calculated only from observability matrix

    [Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA
    y = y+1;
    iv(y) = i;
    OMAfreq = Result.Parameters.NaFreq;
    for j = 1:length(Result.Parameters.NaFreq)
        Acc(j,y) = min(OMAfreq(j),omegas(j))/max(OMAfreq(j),omegas(j));
    end
end
figure
hold on
plot(iv,Acc(1,:))
plot(iv,Acc(2,:))
plot(iv,Acc(3,:))
plot(iv,Acc(4,:))
plot(iv,Acc(5,:))
hold off
title('Convergence analysis of ERA based on natural frequencies')
xlabel('Number of rows in Hankel matrix')
ylabel('Accuracy')
legend('f1','f2','f3','f4','f5')