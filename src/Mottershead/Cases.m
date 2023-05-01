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

% a simple version
Wtheta = eye(3);

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm - K*um;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r

%%
%%%%%%%%%%%%%%
%%% CASE 2 %%%
%%%%%%%%%%%%%%
clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the first measurement
% force vector
fm1 = [1, 0]';

% measured static deflection
um1 = [0.60811;0.27027];

% The sensitivity matrix:
G = [um1(1),um1(1)-um1(2),0;0,-um1(1)+um1(2),um1(2)];

% condition number of G
con1 = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% The parameter weighting matrix for regularization problem from equation
% (19) and (20)
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm1 - K*um1;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
dx1 = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the second measurement
% force vector
fm2 = [0, 1]';

% measured static deflection
um2 = [0.27027;0.67568];

% The sensitivity matrix:
G = [um2(1),um2(1)-um2(2),0;0,-um2(1)+um2(2),um2(2)];

% condition number of G
con2 = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% The parameter weighting matrix for regularization problem from equation
% (19) and (20) 
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% stiffness parameters in the initial model
k1 = 1+dx1(1);
k2 = 1+dx1(2);
k3 = 1+dx1(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm2 - K*um2;  

% the difference with regularization
dx2 = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

dx = dx1 + dx2


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE 2 ill-condition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the first measurement
% force vector
fm1 = [1, 0]';

% measured static deflection
um1 = [0.60811;0.27027]+[-0.0002;0.0003];

% The sensitivity matrix:
G = [um1(1),um1(1)-um1(2),0;0,-um1(1)+um1(2),um1(2)];

% condition number of G
con1 = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% The parameter weighting matrix for regularization problem from equation
% (19) and (20)
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm1 - K*um1;

% plotting to find the regularisation parameter:
lambda = 0.0179;  

% the difference with regularization
dx1 = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the second measurement
% force vector
fm3 = [1.05, 0]';

% measured static deflection
um3 = [0.63851;0.28378]+[0.0010;-0.0002];

% The sensitivity matrix:
G = [um3(1),um3(1)-um3(2),0;0,-um3(1)+um3(2),um3(2)];

% condition number of G
con2 = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% The parameter weighting matrix for regularization problem from equation
% (19) and (20) 
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% stiffness parameters in the initial model
k1 = 1+dx1(1);
k2 = 1+dx1(2);
k3 = 1+dx1(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm3 - K*um3;  

% the difference with regularization
dx2 = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

dx = dx1 + dx2;

k1f = 1+dx(1);
k2f = 1+dx(2);
k3f = 1+dx(3);
kf = [k1f;k2f;k3f]