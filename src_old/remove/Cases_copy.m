clear;clc;close all
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
Weps = [1/2,0;0,7.65/2];

% The parameter weighting matrix for regularization problem
x = diag(G'*Weps*G);
Gamma = [x(1),0,0;0,x(2),0;0,0,x(3)];
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% plotting to find the regularisation parameter:
lambda = 10^-2;  

% the difference with regularization
Deltatheta = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r

%%
%%%%%%%%%%%%%%
%%% JAn %%%
%%%%%%%%%%%%%%

fm = [1 ; 0];


k1m = 1.2;

k2m = 0.8;

k3m = 1;

 

Km = [k1m+k2m -k2m ; -k2m k2m+k3m];

 

um = Km\fm;

 

G = [ um(1) um(1)-um(2) 0 ; 0 -um(1)+um(2) um(2) ];

 

x = [1 ; 1 ; 1];

x0 = x;

 

% weighting matrices

We = diag(3);

Gam = diag(diag(G'*We*G));

iGam = inv(Gam);

Wt = mean(diag(Gam))/mean(diag(Gam))*iGam;

Wt = eye(3);

 

% regularization

lam = 0.1;

 

% for loop

for ii = 1:20

 

    k1 = x(1);

    k2 = x(2);

    k3 = x(3);

 

    K = [k1+k2 -k2 ; -k2 k2+k3];

 

    f = K*um;

 

    dz = fm - f;

 

    SS = G'*We*G + lam^2*Wt;

    RR = G'*We*dz;

    dx = SS\RR;

 

    x = x + dx;

 

end

 

% change in optimization variables

dx = x - x0
