%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame eigenvalue and mode shape residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath(genpath('npy-matlab-master'),genpath('..\..\data'),genpath('..\Costfunctions\functions'))
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);

if q == 1
    % measuered natural frequencies from OMA SSI-cov [Hz]
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_5_2_1.npy');
    
    % Natural frequencies [rad/s]
    SSIomega = SSIFreq * 2 * pi;
    
    % squared
    omegasq_m = SSIomega.^2;
    
    % Load mode shapes
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_5_2_1.npy');
elseif q == 2
    % measuered natural frequencies from OMA SSI-cov [Hz]
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_no_damp.npy');
    
    % Natural frequencies [rad/s]
    SSIomega = SSIFreq * 2 * pi;
    
    % squared
    omegasq_m = SSIomega.^2;
    
    % Load mode shapes
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_no_damp.npy');
end

%  Analytical natural frequencies
% loading the model parameters
filename = load('..\..\data\modelprop.mat'); 

% stiffness parameters
k = filename.k';

% Mass matrix
M = filename.M;

% reshape mode shapes
SSIphi_re = reshape(SSIphi,25,1);

% symmetric weighting matrix
Weps = diag([omegasq_m./max(omegasq_m);SSIphi_re])^-2;

%Weps = eye(30);

% Difine lambda value:
if q == 1
    lambda =  sqrt(0.0984); 
elseif q == 2
    lambda =  sqrt(0.0855); 
end

dx = zeros(length(k),1);

% number of iterations
ni = 100;
% Iterazation over the stiffness
for ii = 1:ni
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
    Umac(:,:,ii) = U;

    % Residual mode shapes
    rmod(:,ii) = SSIphi_re - reshape(U,25,1);

    % Residual eigenvalue
    reig(:,ii) = (omegasq_m - omega_squared);
       
    % Combines residuals
    r(:,ii) = [reig(:,ii);rmod(:,ii)];

    disp(sum(abs(r(:,ii))))

    % New sensitivity matrix
    syms k1 k2 k3 k4 k5 m1 m2 m3 m4 m5
    
    Msym = [m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];
    
    Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];
    
    Geig = zeros(length(omegas),length(k));
    
    % number of modes
    nm = length(omegas);
    % number of parameters
    np = length(k);
    
    pars = [k1,k2,k3,k4,k5];
    parm = [m1,m2,m3,m4,m5];
    for j = 1:nm
        for kk = 1:np
            par = pars(kk);
            Geig(j,kk) = (U(:,j)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j));
        end
    end
    
    % The sensitivity matrix out from the analytical mode shapes:
    syms k1 k2 k3 k4 k5 
    
    Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];
    
    
    % number of modes
    nm = length(SSIomega);
    
    % The sensitivity matrix out from the analytical mode shapes:
    % the a factor
    for j = 1:nm % looping over the mode shape
        for kk = 1:np % looping over the parameters
            for h = 1:nm % Looping over modes shape included in approximated
                 if j == h 
                    par = pars(kk);
                    as(j,kk,h) = -1/2 * U(:,j)'*(diff(Msym,par))*U(:,j);
                elseif j ~= h 
                    par = pars(kk);
                    as(j,kk,h) = (U(:,h)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j))/(omega_squared(j)-omega_squared(h));
                end
               
            end
        end
    end
    
    Gmod = zeros(25,5);
    % now finding the sensivity matrix for the stiffness parameters
    ind = 1;
    b1 = zeros(5,1);
    b = zeros(5,1);
    for j = 1:nm
        for kk = 1:np
            for h = 1:nm
              b1 = as(j,kk,h)*U(:,h);
              b = b+b1;
              if h == nm
                Gmod(ind,kk) = b(1);
                Gmod(ind+1,kk) = b(2);
                Gmod(ind+2,kk) = b(3);
                Gmod(ind+3,kk) = b(4);
                Gmod(ind+4,kk) = b(5);
              end
            end
        end
        ind = ind+5;
    end
    
    % The big sensivity matrix
    G = [Geig;Gmod];
    
    % condition number
    con(ii) = cond(G);

    % a simple version of the parameter weighting matrix
    %Wtheta = eye(size(G,2));
    % The parameter weighting matrix for regularization problem from equation
    % (19) and (20)
    Gamma = diag(diag(G'*Weps*G));
    Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

    
    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r(:,ii);
    GNEW(:,:,ii) = G;
end
% requirement
req = rank(G'*Weps*G) == size(G,2);
disp(req);

% Convergence plot for the relative failure for frequencies
figure (1)
hold on
for j = 1:5
    err(j,:) = abs(r(j,:))/omegasq_m(j)*100;
    plot(err(j,:))
end
legend('Frequency 1', 'Frequency 2','Frequency 3', 'Frequency 4','Frequency 5')
xlabel('Iterations [-]')
ylabel('Relative error [%]')
hold off


% Display accuracy of mode shapes
% CrossMAC plot of mode shapes
for i = 1:ni
    mac=crossMACnm(SSIphi,Umac(:,:,i));
    for j = 1:5
        dmac = diag(mac);
        acc(j,i) = dmac(j)*100;
    end
end
% Convergence plot for the relative failure for mode shapes
figure (2)
hold on
for j = 1:5
    plot(acc(j,:))
end
legend('Mode shape 1', 'Mode shape 2','Mode shape 3', 'Mode shape 4','Mode shape 5')
xlabel('Iterations [-]')
ylabel('MAC [%]')
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

% Save updated system matrices
if q == 1
    save('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_high.mat','knew','Knew','U','fn');
elseif q == 2
    save('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat','knew','Knew','U','fn');
end
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
figure (4)
loglog(Jeps,Jthe)
grid on
% xlim([10^-4 10^-1+10^-1/2])
% ylim([10^-3 10^-0+10^-0/1.5])
xlabel('norm (Residual)')
ylabel('norm (Stiffness Change)')
title('L-curve')
% Finding the optimal value for lambda
Val = 823.941;
index = find(Jeps >= Val,1);
lamopt = lambda(index);
