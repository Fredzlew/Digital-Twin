clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%Model Parameters and excitation
%--------------------------------------------------------------------------
% Set global random seed
rng(2)

filename = load('modelprop.mat'); % Loads mass and stiffness matrices
omega_min = 1.70; % Minimum natural frequency
zeta_min = 0.015; % Minimum threshold for desired damping ratio
alpha = zeta_min*omega_min; % Rayleigh damping coefficient
beta = zeta_min/omega_min; % Rayleigh damping coefficient
% Known system matrices
M=filename.M; % Mass matrix
K=filename.K; % Stiffness matrix
C=0;%alpha*M+beta*K;
f=5*randn(5,400000);
fs=100;

%Apply modal superposition to get response
%--------------------------------------------------------------------------

n=size(f,1);
dt=1/fs; %sampling rate
[Us, Values]=eig(K,M);
Freq=sqrt(diag(Values))/(2*pi); % undamped natural frequency
steps=size(f,2);

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

Mn=diag(Vectors'*M*Vectors); % uncoupled mass
Cn=diag(Vectors'*C*Vectors); % uncoupled damping
Kn=diag(Vectors'*K*Vectors); % uncoupled stifness
wn=sqrt(diag(Values));
zeta=Cn./(sqrt(2.*Mn.*Kn));  % damping ratio
wd=wn.*sqrt(1-zeta.^2);

fn=Vectors'*f; % generalized input force matrix

t=[0:dt:dt*steps-dt];

for i=1:1:n
    
    h(i,:)=(1/(Mn(i)*wd(i))).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t); %transfer function of displacement
    hd(i,:)=(1/(Mn(i)*wd(i))).*(-zeta(i).*wn(i).*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)+wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)); %transfer function of velocity
    hdd(i,:)=(1/(Mn(i)*wd(i))).*((zeta(i).*wn(i))^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)-zeta(i).*wn(i).*wd(i).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t)-wd(i).*((zeta(i).*wn(i)).*exp(-zeta(i)*wn(i)*t).*cos(wd(i)*t))-wd(i)^2.*exp(-zeta(i)*wn(i)*t).*sin(wd(i)*t)); %transfer function of acceleration
    
    qq=conv(fn(i,:),h(i,:))*dt;
    qqd=conv(fn(i,:),hd(i,:))*dt;
    qqdd=conv(fn(i,:),hdd(i,:))*dt;
    
    q(i,:)=qq(1:steps); % modal displacement
    qd(i,:)=qqd(1:steps); % modal velocity
    qdd(i,:)=qqdd(1:steps); % modal acceleration
       
end

x=Vectors*q; %displacement
v=Vectors*qd; %vecloity
a=Vectors*qdd; %vecloity

%Add noise to excitation and response
%--------------------------------------------------------------------------
rng(3)
f2=f+0.05*randn(5,length(f));
a2=a+0.05*randn(5,length(f));
v2=v+0.05*randn(5,length(f));
dis=x+0.05*randn(5,length(f));

save('.\data\data_sim.mat','dis');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% NEWMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions
x0 = zeros(size(f,1),1);
v0 = x0;

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
dis_new=x_new+0.05*randn(5,length(f));

save('.\data\data_sim_newmark.mat','dis_new');