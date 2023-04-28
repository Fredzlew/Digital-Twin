%%%%%%%%%%%%%%
%%% CASE 1 %%%
%%%%%%%%%%%%%%

% force vector
fm = [1, 0]';

% measured static deflection
um = [0.60811;0.27027];

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm - K*um;

% The sensitivity matrix:
G = [um(1),um(1)-um(2),0;0,-um(1)+um(2),um(2)];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = (diag(fm))^-2;
Weps = [1,0;0,7.65]
Weps = [1/2,0;0,7.65/2]
%Weps = (1/norm(r)^2)*(r*r') 

% The parameter weighting matrix for regularization problem
x = diag(G'*Weps*G);
Gamma = [x(1),0,0;0,x(2),0;0,0,x(3)];
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
Deltatheta = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r

% The difference with parameter estimation
%Deltatheta = (G'*Weps*G)^-1*G'*Weps*r;

%%
%%%%%%%%%%%%%%
%%% CASE 2 %%%
%%%%%%%%%%%%%%

% force vector
fm = [0, 1]';

% measured static deflection
um = [0.27027;0.67568];

% stiffness parameters in the initial model
k1 = 1+Deltatheta(1);
k2 = 1+Deltatheta(2);
k3 = 1+Deltatheta(3);

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm - K*um;

% The sensitivity matrix:
G = [um(1),um(1)-um(2),0;0,-um(1)+um(2),um(2)];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = [1/2,0;0,7.65/2]
%Weps = (1/norm(r)^2)*(r*r') 

% The parameter weighting matrix for regularization problem
x = diag(G'*Weps*G);
Gamma = [x(1),0,0;0,x(2),0;0,0,x(3)];
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
Deltatheta = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r

% The difference with parameter estimation
%Deltatheta = (G'*Weps*G)^-1*G'*Weps*r;
%%
%%%%%%%%%%%%%%
%%% CASE 2 %%%
%%%%%%%%%%%%%%
% force vector
fm1 = [1, 0]';
fm2 = [0, 1]';
fm = [fm1,fm2];


% measured static deflection
um1 = [0.60811;0.27027];
um2 = [0.27027;0.67568];
um = [um1,um2];

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% input force residual
r = fm - K*um;

% The sensitivity matrix:
G = [um(1),um(1)-um(2),0;0,-um(1)+um(2),um(2)];



% condition number of G
con = cond(G);

% symmetric weighting matrix
%Weps = (diag(fm))^-2;
Weps = [1,0;0,1];
%Weps = (1/norm(r)^2)*(r*r') 

% The parameter weighting matrix for regularization problem
x = diag(G'*Weps*G);
Gamma = [x(1),0,0;0,x(2),0;0,0,x(3)];
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
Deltatheta = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

% The difference with parameter estimation
Deltatheta = (G'*Weps*G)^-1*G'*Weps*r