%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame eigenvalue residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath(genpath('npy-matlab-master'),genpath('..\..\data'))
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

% symmetric weighting matrix
Weps = diag(omegasq_m./max(omegasq_m))^-2;
%Weps = eye(size(G,1));


% Difine lambda value:
if q == 1
    lambda =  sqrt(0.3385); 
elseif q == 2
    lambda =  sqrt(0.2779); 
end

% initial conditions
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
    
    % natural frequencies [rad/s]
    omegas = omega(iw);
    
    % natural frequencies squared
    omega_squared = omegas.^2.;

    % Residual
    r(:,ii) = (omegasq_m - omega_squared);


    % CALCULATE A NEW G
    syms k1 k2 k3 k4 k5 m1 m2 m3 m4 m5
    
    Msym = [m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];
    
    Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];
    
    G = zeros(length(omegas),length(k));
    
    % number of modes
    nm = length(omegas);
    % number of parameters
    np = length(k);
    
    pars = [k1,k2,k3,k4,k5];
    parm = [m1,m2,m3,m4,m5];
    for j = 1:nm
        for kk = 1:np
            par = pars(kk);
            G(j,kk) = (U(:,j)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j));
        end
    end
    GNEW(:,:,ii) = G;

    % condition number of G
    con(ii) = cond(G);

    disp(sum(abs(r(:,ii)))) 

     
    % The parameter weighting matrix for regularization problem from equation
    % (19) and (20)
    Gamma = diag(diag(G'*Weps*G));
    Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;
    % Wtheta = eye(size(G,2));
    
    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r(:,ii);
end


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
    figure (1)
    loglog(Jeps,Jthe)
    grid on
    xlabel('norm (Residual)','FontSize',14)
    ylabel('norm (Stiffness Change)','FontSize',14)
    title('L-curve','FontSize',20)
    % Finding the optimal value for lambda
    Val = 154.603; % X-vlaue
    index = find(Jeps >= Val,1);
    lamopt = lambda(index);
    if q == 1
        save('.\data_updated_par_sens\Eigenvalue_residual_L_curve_high.mat','Jeps','Jthe');
    elseif q == 2
        save('.\data_updated_par_sens\Eigenvalue_residual_L_curve_no_damp.mat','Jeps','Jthe');
    end
end

if ni == 100
    if q == 1
        save('.\data_updated_par_sens\Eigenvalue_residual_high.mat','knew','Knew','U','fn','err');
    elseif q == 2
        save('.\data_updated_par_sens\Eigenvalue_residual_no_damp.mat','knew','Knew','U','fn','err');
    end
end
