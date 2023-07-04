%%%%%%%%%%%%%%
%%% CASE 1 %%%
%%%%%%%%%%%%%%
clear all
clc
% force vector
omega_m = [1.1, 1.8]';

% Mass matrix
M = [1,0;0,1];

% stiffness parameters in the initial model
k1 = 1;
k2 = 1;
k3 = 1;

% stiffness matrix
K = [k1+k2,-k2;-k2,k2+k3];

% eigenvalue problem
[U,D] = eig(K,M);

% natural frequencies from eigenvalues
omega = real(sqrt(diag(D)));

% sort frequencies and mode shapes
[~,iw] = sort(omega);

% natural frequencies [rad/s]
omegas = omega(iw);

% mode shapes
Us = U(:,iw);

% normalization
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:length(omegas)
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:length(omegas)
        U(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization


% The sensitivity matrix:
G = [U(1,1)^2,(U(1,1)-U(2,1))^2,U(2,1)^2;U(1,2)^2,(U(1,2)-U(2,2))^2,U(2,2)^2];

% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(2);

% a simple version of the parameter weighting matrix
Wtheta = eye(3);

% input force residual
r = omega_m - omegas;

% Difne lambda value:
lambda = 10^-2;  

% the difference with regularization
dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r


