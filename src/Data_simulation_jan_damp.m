clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Set global random seed
rng(2)

filename = load('modelprop_jan.mat'); % Loads mass and stiffness matrices
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
[Us, Values]=eig(K,M); % Solve eigenvalue problem
omegas=sqrt(diag(Values)); % undamped natural frequency
omega_min = omegas(1); % Minimum natural frequency [rad/s]
zeta_min = 0.002; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
betas = zeta_min/omega_min; % Rayleigh damping coefficient
C=alpha*M+betas*K; % Damping matrix
zetas = (alpha./omegas + betas.*omegas)./2; % Calculate damping ratios
f=2e1*randn(5,2.3e6);  
 
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
dis_new=x_new+0.03*x_new.*randi([-1,1],size(x_new,1),size(x_new,2));
save('.\data\data_sim_newmark_jan_damp.mat','dis_new','zetas');
save('.\python\data\data_sim_newmark_jan_damp.mat','dis_new','zetas');

% Calculates the periodic error for Average Acceleration
disp(1/12*(omegas(5)*dt)^2)
