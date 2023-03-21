clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Set global random seed
rng(2)

filename = load('modelprop_jan.mat'); % Loads mass and stiffness matrices
omega_min = 1.7474; % Minimum natural frequency
zeta_min = 0.002; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient
% Known system matrices
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=alpha*M+beta*K; % Damping matrix
zeta=C./(2*sqrt(M.*K)); % Calculate damping ratios
f=1e2*randn(5,2.3e6);  
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% NEWMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions
x0 = zeros(size(f,1),1);
v0 = x0;

% Sampling rate
dt = 0.001;

% Number of time steps
N = length(f)-1;

% Average Acceleration
gamma = 0.5;
beta = 0.25;
% Fox-Goodwin
% gamma = 0.5;
% beta = 1/12;

% Newmark time integration
[x_new,v_new,a_new,t_new] = newmark(K,C,M,x0,v0,...
    dt,N,f,beta,gamma);

% Add noise to response
dis_new=x_new+0.005*randn(5,length(f));
save('.\data\data_sim_newmark_jan_damp.mat','dis_new','zeta');
save('.\python\data\data_sim_newmark_jan_damp.mat','dis_new','zeta');
