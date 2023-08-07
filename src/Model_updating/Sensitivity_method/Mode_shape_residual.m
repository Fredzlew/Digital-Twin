%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame mode_shape residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath(genpath('..\costfunctions\functions'),genpath('..\..\npy-matlab-master'),genpath('..\..\data'))
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);

if q == 1
    % measuered natural frequencies from OMA SSI-cov [Hz]
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_5_2_1.npy');
    
    % Load mode shapes
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_5_2_1.npy');
elseif q == 2
    % measuered natural frequencies from OMA SSI-cov [Hz]
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_no_damp.npy');
    
    % Load mode shapes
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par\SSImodes_no_damp.npy');
end

%  Analytical natural frequencies
% loading the model parameters
filename = load('..\..\data\modelprop.mat'); 

% load stiffness parameters
k = filename.k';

% Mass matrix
M = filename.M;

% reshape mode shapes
SSIphi_re = reshape(SSIphi,25,1);

% symmetric weighting matrix
%Weps = eye(size(G,1));
Weps = diag(SSIphi_re)^-2;

% Difine lambda value:
if q == 1
    lambda =  sqrt(1.0001);  
elseif q == 2
    lambda =  sqrt(1.4001); 
end

% initial condition
dx = zeros(length(k),1);

% number of iterations
ni = 1;

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
    for j = 1:length(omega)
        if abs(MVec_x(j)) > abs(mVec_x(j))
            mxVec_x(j) = MVec_x(j);
        else
            mxVec_x(j) = mVec_x(j);
        end
        for l = 1:length(omega)
            U(l,j) = Us(l,j)/mxVec_x(j);
        end
    end % end normalization
    Umac(:,:,ii) = U;
    % Residual
    r(:,ii) = (SSIphi_re) - (reshape(U,25,1));
    disp(sum(abs(r(:,ii))))

    % New sensitivity matrix
    syms k1 k2 k3 k4 k5 m1 m2 m3 m4 m5
    
    Msym = [m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];
    
    Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];
    
    
    % number of modes
    nm = length(omega);
    
    % number of parameters
    np = length(K);
    
    pars = [k1,k2,k3,k4,k5];
    parm = [m1,m2,m3,m4,m5];
    
    % The sensitivity matrix out from the analytical mode shapes:
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
    
    G = zeros(25,5);
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
                G(ind,kk) = b(1);
                G(ind+1,kk) = b(2);
                G(ind+2,kk) = b(3);
                G(ind+3,kk) = b(4);
                G(ind+4,kk) = b(5);
              end
            end
        end
        ind = ind+5;
    end
    % condition number of G
    con(ii) = cond(G);

    % a simple version of the parameter weighting matrix
    %Wtheta = eye(size(G,2));
    % The parameter weighting matrix for regularization problem from equation
    % (19) and (20)
    Gamma = diag(diag(G'*Weps*G));
    Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r(:,ii);

end

% Convergence plot

% Display accuracy of mode shapes
% CrossMAC plot of mode shapes
for i = 1:ni
    mac=crossMAC(SSIphi,Umac(:,:,i));
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
fn = omegasnew/(2*pi);

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


%% Plotting L curve only for the first iteration

% plotting
if ni == 1
    qq = 1;
    for i = linspace(0.00000001,1000,10000)
        lambda(qq) = i;
        dx = ((G'*Weps*G)+(lambda(qq)^2*Wtheta))^(-1)*G'*Weps*r(:,end);
        eps = r(:,end)-G*dx; 
        Jeps(qq) = sqrt(eps'*Weps*eps);
        Jthe(qq)  = sqrt(dx'*Wtheta*dx);
        qq = qq + 1;
    end
    % plotting the L curve
    figure (4)
    loglog(Jeps,Jthe)
    grid on
    xlabel('norm (Residual)')
    ylabel('norm (Stiffness Change)')
    title('L-curve')
    % Finding the optimal value for lambda
    Val = 0.748157;
    index = find(Jeps >= Val,1);
    lamopt = lambda(index);
    if q == 1
        save('.\data_updated_par_sens\Mode_shape_residual_L_curve_high.mat','Jeps','Jthe');
    elseif q == 2
        save('.\data_updated_par_sens\Mode_shape_residual_L_curve_no_damp.mat','Jeps','Jthe');
    end
end

if ni == 100
    if q == 1
        save('.\data_updated_par_sens\Mode_shape_residual_high.mat','knew','Knew','U','fn','acc','con');
    elseif q == 2
        save('.\data_updated_par_sens\Mode_shape_residual_no_damp.mat','knew','Knew','U','fn','acc','con');
    end
end