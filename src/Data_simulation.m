clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
%for i = 1:1000
% Set global random seed
rng(1)
i = 1;
filename = load('.\data\modelprop_jan.mat'); % Loads mass and stiffness matrices
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
[Us, Values]=eig(K,M); % Solve eigenvalue problem
omegas=sqrt(diag(Values)); % undamped natural frequency
% normalization
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        U(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization
C=0; % Damping matrix
f=2e1*randn(1,1.2e7).*ones(5,1);  
 
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
dis_new=x_new; %+0.03*x_new.*randi([-1,1],size(x_new,1),size(x_new,2));

%save('.\data\data_sim_newmark_jan_damp.mat','dis_new','zetas');
%filename2 = sprintf('python\data\%04d_data_sim_newmark_jan_damp.mat',i);
save(['data\' num2str(i) '_data_sim_newmark_jan.mat'],'dis_new','U');
% Periodic error
PERFEJL = 1/12*(omegas(5)*dt)^2;

%end