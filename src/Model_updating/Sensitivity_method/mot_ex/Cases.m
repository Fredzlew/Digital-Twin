%%%%%%%%%%%%%%
%%% CASE 1 %%%
%%%%%%%%%%%%%%
clear all
clc
% force vector
fm = [1, 0]';

% measured static deflection
um = [0.60811;0.27027];

% The sensitivity matrix:
G = [um(1),um(1)-um(2),0;0,-um(1)+um(2),um(2)];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% a simple version of the parameter weighting matrix
Wtheta = eye(3);

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

dx = zeros(3,1);
dxt = dx;
k0 = [k1;k2;k3];

for ii=1:1
% The updated stiffnesses
k1 = k1+dx(1);
k2 = k2+dx(2);
k3 = k3+dx(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm - K*um;

% Define lambda value:
lambda = 10^-2;  

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

% Total parameter change
dxt = dxt + dx;
end

disp(dxt)
kf = k0+dxt

%%
%%%%%%%%%%%%%%
%%% CASE 2 %%%
%%%%%%%%%%%%%%
clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% force vectors
fm1 = [1, 0]';
fm2 = [0, 1]';
fm = [fm1;fm2];

% measured static deflections
um1 = [0.60811;0.27027];
um2 = [0.27027;0.67568];

% The sensitivity matrices:
G1 = [um1(1),um1(1)-um1(2),0;0,-um1(1)+um1(2),um1(2)];
G2 = [um2(1),um2(1)-um2(2),0;0,-um2(1)+um2(2),um2(2)];
G = [G1;G2];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(size(G,1));

% a simple version of the parameter weighting matrix
Wtheta = eye(size(G,2));

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

dx = zeros(3,1);
dxt = dx;
k0 = [k1;k2;k3];

for ii=1:1
% The updated stiffnesses
k1 = k1+dx(1);
k2 = k2+dx(2);
k3 = k3+dx(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% Analytical forces
f1 = K*um1;
f2 = K*um2;
f = [f1;f2];

% input force residual
r = fm - f;

% plotting to find the regularisation parameter:
lambda = 10^-3; 

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

% Total parameter change
dxt = dxt + dx;
end

disp(dxt)
kf = k0+dxt

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 2 ill-condition (Identity matrix) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% force vectors
fm1 = [1, 0]';
fm3 = [1.05, 0]';
fm = [fm1;fm3];

% measured static deflections
um1 = [0.60811;0.27027]+[-0.0002;0.0003];
um3 = [0.63851;0.28378]+[0.0010;-0.0002];

% The sensitivity matrices:
G1 = [um1(1),um1(1)-um1(2),0;0,-um1(1)+um1(2),um1(2)];
G3 = [um3(1),um3(1)-um3(2),0;0,-um3(1)+um3(2),um3(2)];
G = [G1;G3];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(size(G,1));

% a simple version of the parameter weighting matrix
Wtheta = eye(size(G,2));

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

dx = zeros(3,1);
dxt = dx;
k0 = [k1;k2;k3];

for ii = 1:1
% The updated stiffnesses
k1 = k1+dx(1);
k2 = k2+dx(2);
k3 = k3+dx(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% Analytical forces
f1 = K*um1;
f3 = K*um3;
f = [f1;f3];

% input force residual
r = fm - f;

% Define lambda value
lambda = 0.0093; 

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

% Total parameter change
dxt = dxt + dx;
end

disp(dxt)
kf = k0+dxt

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 2 ill-condition (eq. 19 and 20) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% force vectors
fm1 = [1, 0]';
fm3 = [1.05, 0]';
fm = [fm1;fm3];

% measured static deflections
um1 = [0.60811;0.27027]+[-0.0002;0.0003];
um3 = [0.63851;0.28378]+[0.0010;-0.0002];

% The sensitivity matrices:
G1 = [um1(1),um1(1)-um1(2),0;0,-um1(1)+um1(2),um1(2)];
G3 = [um3(1),um3(1)-um3(2),0;0,-um3(1)+um3(2),um3(2)];
G = [G1;G3];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(size(G,1));

% The parameter weighting matrix for regularization problem from equation
% (19) and (20)
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

dx = zeros(3,1);
dxt = dx;
k0 = [k1;k2;k3];

for ii = 1:1
% The updated stiffnesses
k1 = k1+dx(1);
k2 = k2+dx(2);
k3 = k3+dx(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% Analytical forces
f1 = K*um1;
f3 = K*um3;
f = [f1;f3];

% input force residual
r = fm - f;

% Define lambda value
lambda = 0.0179; 

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

% Total parameter change
dxt = dxt + dx;
end

disp(dxt)
kf = k0+dxt

%% Plotting L curve

% plotting
q = 1;
for i = linspace(0.0001,10,100000)
    lambda = i;
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;
    eps = r-G*dx;
    Jeps(q) = sqrt(eps'*Weps*eps);
    Jthe(q)  = sqrt(dx'*Wtheta*dx);
    q = q + 1;
end
% plotting the L curve
figure (1)
loglog(Jeps,Jthe)
grid on
%xlim([10^-4 10^-1+10^-1/2])
%ylim([10^-3 10^-0+10^-0/1.5])
xlabel('norm (Residual)')
ylabel('norm (Stiffness Change)')
title('L-curve')

% plotting the  norm to the regularization parameter
% lambda square:
%lam2 = linspace(10^-10,10^0,100000);


%q = 1;
%for ii = lam2
    %stiffdx(:,q) = ((G'*Weps*G)+(ii*Wtheta))^(-1)*G'*Weps*r;
    %q = q + 1;
%end

% plotting the  stiffness change to the reqularization parameter
%figure (2)
%semilogx(lam2,stiffdx(1,:),lam2,stiffdx(2,:),lam2,stiffdx(3,:))
%xlim([10^-10 10^0])
%ylim([-1.5 1.5])
%grid on
%xlabel('Regularization Parameter, lambda^2')
%ylabel('Stiffness Change')