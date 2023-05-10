%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame mode_shape residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('functions'),genpath('OMA'),genpath('python'),genpath('npy-matlab-master'),genpath('data'),genpath('Modal_parameters_anela'))

% measuered natural frequencies from OMA SSI-cov 
% SSIFreq = readNPY('SSIomega_5_2_1.npy');
SSIFreq =  [1.6570; 5.0168; 7.8984; 10.1144; 11.5872];

SSIomega = SSIFreq * 2 * pi;
% squared
omegasq_m = SSIomega.^2;

%SSIphi = readNPY('SSIphi_5_2_1.npy');
SSIphi = [0.3164, 0.7748, 1.0000, 1.0000, -0.4440;...
    0.6301, 1.0000, 0.2081, -0.9371,0.8281;...
    0.7783, 0.4530, -0.7971, -0.0186, -0.9820;...
    1.0000, -0.3185, -0.3893, 0.8750, 1.0000;...
    0.9923, -0.7864, 0.6152, -0.5075, -0.3861];



%  Analytical natural frequencies
% loading the model parameters
% filename = load('modelprop_jan.mat'); 

% stiffness parameters
% k = filename.k';
k = [4.0700; 3.2289; 3.3756; 3.5224; 3.6740]*1e3;

% stiffness matrix (stiffness parameters in the initial model)
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);

% Mass matrix
% M = filename.M;

M = [2.3553, 0, 0, 0, 0;...
     0, 2.3690, 0, 0, 0;...
     0, 0, 2.3690, 0, 0;...
     0, 0, 0, 2.3690, 0;...
     0, 0, 0, 0, 2.4467];

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

% mode shapes for the analytical model
Us = U(:,iw);

% normalization of the mode shapes
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
% syms k1 k2 k3 k4 k5 

% K = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];
% 
% % number of modes
% nm = length(SSIomega);
% for j = 1:nm
%     for kk = 1:nm
%         for h = 1:nm
%             if j == h
%                 a(j,kk,h) = 0;
%             elseif kk == 1 
%                 stif = k1;
%                 a(j,kk,h) = (SSIphi(:,h)'*(diff(K,stif))*SSIphi(:,j))/(omegasq_m(j)-omegasq_m(h));
%                 b(:,kk,h) = a(j,kk,h)*SSIphi(:,h);
%             elseif kk == 2
%                 stif = k2;
%                 a(j,kk,h) = (SSIphi(:,h)'*(diff(K,stif))*SSIphi(:,j))/(omegasq_m(j)-omegasq_m(h));
%                 b(:,kk,h) = a(j,kk,h)*SSIphi(:,h);
%             elseif kk == 3
%                 stif = k3;
%                 a(j,kk,h) = (SSIphi(:,h)'*(diff(K,stif))*SSIphi(:,j))/(omegasq_m(j)-omegasq_m(h));
%                 b(:,kk,h) = a(j,kk,h)*SSIphi(:,h);
%             elseif kk == 4
%                 stif = k4;
%                 a(j,kk,h) = (SSIphi(:,h)'*(diff(K,stif))*SSIphi(:,j))/(omegasq_m(j)-omegasq_m(h));
%                 b(:,kk,h) = a(j,kk,h)*SSIphi(:,h);
%             elseif kk == 5
%                 stif = k5;
%                 a(j,kk,h) = (SSIphi(:,h)'*(diff(K,stif))*SSIphi(:,j))/(omegasq_m(j)-omegasq_m(h));
%                 b(:,kk,h) = a(j,kk,h)*SSIphi(:,h);
%             end
%            
%         end
%     end
% end
% for j = 1:nm
%     for kk = 1:nm
%         G(j,kk) = sum(b(j,kk,:));
%     end
% end

% The sensitivity matrix out from the analytical mode shapes:
syms k1 k2 k3 k4 k5 

Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];


% number of modes
nm = length(SSIomega);

% making the a factor
for j = 1:nm
    for kk = 1:nm
        for h = 1:nm
            if j == h
                a(j,kk,h) = 0;
            elseif kk == 1 
                stif = k1;
                a(j,kk,h) = (U(:,h)'*(diff(Ksym,stif))*U(:,j))/(omega_squared(j)-omega_squared(h));
            elseif kk == 2
                stif = k2;
                a(j,kk,h) = (U(:,h)'*(diff(Ksym,stif))*U(:,j))/(omega_squared(j)-omega_squared(h));
            elseif kk == 3
                stif = k3;
                a(j,kk,h) = (U(:,h)'*(diff(Ksym,stif))*U(:,j))/(omega_squared(j)-omega_squared(h));
            elseif kk == 4
                stif = k4;
                a(j,kk,h) = (U(:,h)'*(diff(Ksym,stif))*U(:,j))/(omega_squared(j)-omega_squared(h));
            elseif kk == 5
                stif = k5;
                a(j,kk,h) = (U(:,h)'*(diff(Ksym,stif))*U(:,j))/(omega_squared(j)-omega_squared(h));
            end
           
        end
    end
end

G = zeros(25,5);
% now finding the sensivity matrix
m = 1;
b1 = zeros(5,1);
b = zeros(5,1);
for j = 1:nm
    for kk = 1:nm
        for h = 1:nm
          b1 = a(j,kk,h)*U(:,h);
          b = b+b1;
          if h == nm
            G(m,kk) = b(1);
            G(m+1,kk) = b(2);
            G(m+2,kk) = b(3);
            G(m+3,kk) = b(4);
            G(m+4,kk) = b(5);
          end
        end
    end
    m = m+5;
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

dx = zeros(length(k),1);
SSIphi_re = reshape(SSIphi,25,1);
% Iterazation over the stiffness
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

    % mode shapes for the analytical model
    Us = U(:,iw);
    
    % normalization of the mode shapes
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
    
    % Residual
    r(:,ii) = SSIphi_re - reshape(U,25,1);
    disp(sum(abs(r(:,ii))))

    % Difne lambda value:
    lambda =  sqrt(0.0143);  
    
    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r(:,ii);

end

% Convergence plot

figure
for j = 1:25
    hold on
    plot((r(j,:)))
end
legend('Frequency 1', 'Frequency 2','Frequency 3', 'Frequency 4','Frequency 5')
xlabel('Iterations [-]')
ylabel('Residual')
hold off
%% Calculating new frequencies with the stiffness changes

% stiffness parameters
knew = k+dx;

% Stiffness amtrix
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

%% Plotting L curve only for the first iteration

% plotting
q = 1;
for i = linspace(0.0000000000000000001,100,1000000)
    lambda(q) = i;
    dx = ((G'*Weps*G)+(lambda(q)^2*Wtheta))^(-1)*G'*Weps*r(:,end);
    eps = r(:,end)-G*dx; % fortegn + eller -?
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
% Finding the optimal value for lambda
Val = 6.79461;
index = find(Jeps >= Val,1);
lamopt = lambda(index);
% plotting the  norm to the regularization parameter
% lambda square:
% lam2 = linspace(10^-10,10^0,100000);
% 
% 
% q = 1;
% for ii = lam2
%     stiffdx(:,q) = ((G'*Weps*G)+(ii*Wtheta))^(-1)*G'*Weps*r(:,end);
%     q = q + 1;
% end
% 
% % plotting the  stiffness change to the reqularization parameter
% figure (2)
% semilogx(lam2,stiffdx(1,:),lam2,stiffdx(2,:),lam2,stiffdx(3,:))
% xlim([10^-10 10^0])
% ylim([-1.5 1.5])
% grid on
% xlabel('Regularization Parameter, lambda^2')
% ylabel('Stiffness Change')