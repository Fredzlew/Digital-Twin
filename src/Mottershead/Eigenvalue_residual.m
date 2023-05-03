%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CASE eigenvalue residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
addpath(genpath('functions'),genpath('OMA'),genpath('python'),genpath('npy-matlab-master'),genpath('data'),genpath('Modal_parameters_anela'))

% measuered natural frequencies from OMA SSI-cov and are squared
SSIFreq = readNPY('SSIomega_5_2_1.npy');
SSIomega = SSIFreq * 2 * pi;
SSIphi = readNPY('SSIphi_5_2_1.npy');

omegasq_m = SSIomega.^2;

%  Analytical way to find the natural frequencies
filename = load('modelprop_jan.mat'); % Model parameters from the analytical way

% stiffness parameters
k = filename.k;

% stiffness matrix (stiffness parameters in the initial model)
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);

% Mass matrix
M = filename.M;

% eigenvalue problem
[U,D] = eig(K,M);

% natural frequencies from eigenvalues
omega = real(sqrt(diag(D)));

% sort frequencies and mode shapes
[~,iw] = sort(omega);

% natural frequencies [rad/s]
omegas = omega(iw);

% natural frequencies squared
omega_squared = omegas.^2.;


% mode shapes for the analytical
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

% The sensitivity matrix out from the OMA analysis:
G = zeros(length(omegas),length(k)-1);
for i = 1:length(k)-1
    for j = 1:length(omegas)
        if i == 1
            G(j,i) = SSIphi(j,i)^2;
        else
            G(j,i) = (SSIphi(j,i-1)-SSIphi(j,i))^2;
        end
    end
end
% condition number of G
con = cond(G);

% symmetric weighting matrix
Weps = eye(size(G,1));

% a simple version of the parameter weighting matrix
%Wtheta = eye(size(G,2));
% The parameter weighting matrix for regularization problem from equation
% (19) and (20)
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

dx = zeros(length(k));

for ii = 1:1
    k = k+dx;

    for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
    end
    K(5,5) = k(5);
    
    % eigenvalue problem
    [U,D] = eig(K,M);
    
    % natural frequencies from eigenvalues
    omega = real(sqrt(diag(D)));
    
    % sort frequencies and mode shapes
    [~,iw] = sort(omega);
    
    % natural frequencies [rad/s]
    omegas = omega(iw);
    
    % natural frequencies squared
    omega_squared = omegas.^2.;

    % input force residual
    r = omegasq_m - omega_squared;
    
    % Difne lambda value:
    lambda =  sqrt(0.0719);  
    
    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r;

end

disp(dx)
%%

knew = k(1:4)+dx;

for i = 1:4
    Knew(i,i) = knew(i)+knew(i+1);
    Knew(i,i+1) = -knew(i+1);
    Knew(i+1,i) = -knew(i+1);
end
Knew(5,5) = k(5);

% eigenvalue problem
[U,D] = eig(Knew,M);

% natural frequencies from eigenvalues
omega = real(sqrt(diag(D)));

% sort frequencies and mode shapes
[~,iw] = sort(omega);

% natural frequencies [rad/s]
omegasnew = omega(iw);

% frequencies
fn = omegasnew/(2*pi)

%% Plotting L curve

% plotting
q = 1;
for i = linspace(0.0000000000000000001,100,1000000)
    lambda(q) = i;
    dx = ((G'*Weps*G)+(lambda(q)^2*Wtheta))^(-1)*G'*Weps*r;
    eps = r-G*dx; % fortegn + eller -?
    Jeps(q) = sqrt(eps'*Weps*eps);
    Jthe(q)  = sqrt(dx'*Wtheta*dx);
    q = q + 1;
end
% plotting the L curve
figure (1)
loglog(Jeps,Jthe)
grid on
xlim([10^-4 10^-1+10^-1/2])
ylim([10^-3 10^-0+10^-0/1.5])
xlabel('norm (Residual)')
ylabel('norm (Stiffness Change)')
title('L-curve')

% plotting the  norm to the regularization parameter
% lambda square:
lam2 = linspace(10^-10,10^0,100000);


q = 1;
for ii = lam2
    stiffdx(:,q) = ((G'*Weps*G)+(ii*Wtheta))^(-1)*G'*Weps*r;
    q = q + 1;
end

% plotting the  stiffness change to the reqularization parameter
figure (2)
semilogx(lam2,stiffdx(1,:),lam2,stiffdx(2,:),lam2,stiffdx(3,:))
xlim([10^-10 10^0])
ylim([-1.5 1.5])
grid on
xlabel('Regularization Parameter, lambda^2')
ylabel('Stiffness Change')