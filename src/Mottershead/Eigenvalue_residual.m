%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame eigenvalue residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
addpath(genpath('functions'),genpath('OMA'),genpath('python'),genpath('npy-matlab-master'),genpath('data'),genpath('Modal_parameters_anela'))

% measuered natural frequencies from OMA SSI-cov 
SSIFreq = readNPY('SSIomega_5_2_1.npy');
% SSIFreq =  [1.6570; 5.0168; 7.8984; 10.1144; 11.5872];

SSIomega = SSIFreq * 2 * pi;
% squared
omegasq_m = SSIomega.^2;

SSIphi = readNPY('SSIphi_5_2_1.npy');
% SSIphi = [0.3164, 0.7748, 1.0000, 1.0000, -0.4440;...
%     0.6301, 1.0000, 0.2081, -0.9371,0.8281;...
%     0.7783, 0.4530, -0.7971, -0.0186, -0.9820;...
%     1.0000, -0.3185, -0.3893, 0.8750, 1.0000;...
%     0.9923, -0.7864, 0.6152, -0.5075, -0.3861];



%  Analytical natural frequencies
% loading the model parameters
filename = load('modelprop_jan.mat'); 

% stiffness parameters
k = filename.k';
% k = [4.0700; 3.2289; 3.3756; 3.5224; 3.6740]*1e3;

% stiffness matrix (stiffness parameters in the initial model)
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);

% Mass matrix
M = filename.M;

% M = [2.3553, 0, 0, 0, 0;...
%      0, 2.3690, 0, 0, 0;...
%      0, 0, 2.3690, 0, 0;...
%      0, 0, 0, 2.3690, 0;...
%      0, 0, 0, 0, 2.4467];

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
% G = zeros(length(omegas),length(k));
% for i = 1:length(k)
%     for j = 1:length(omegas)
%         if i == 1
%             G(j,i) = SSIphi(j,i)^2;
%         else
%             G(j,i) = (SSIphi(j,i-1)-SSIphi(j,i))^2;
%         end
%     end
% end

% The sensitivity matrix out from the analytical mode shapes:
G = zeros(length(omegas),length(k));
for i = 1:length(k)
    for j = 1:length(omegas)
        if i == 1
            G(j,i) = U(j,i)^2;
        else
            G(j,i) = (U(j,i-1)-U(j,i))^2;
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

dx = zeros(length(k),1);

% Iterazation over the stiffness
for ii = 1:100
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

    % Residual
    r(:,ii) = omegasq_m - omega_squared;
    disp(sum(abs(r(:,ii))))
    % Difne lambda value:
    lambda =  sqrt(0.0143);  
    
    % the difference with regularization
    dx = ((G'*Weps*G)+(lambda^2*Wtheta))^(-1)*G'*Weps*r(:,ii);

end

% Convergence plot

figure
for j = 1:5
    hold on
    plot((r(j,:)))
end
legend('Frequency 1', 'Frequency 2','Frequency 3', 'Frequency 4','Frequency 5')
xlabel('Iterations [-]','FontSize',14)
ylabel('Residual','FontSize',14)
title('Convergence plot','FontSize',20)
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

% plotting the mode shapes
% dimensions in meters
t = 0.015; % floor height [m]
h = 1*10^-3; % short side of column [m]
b = 30*10^-3; % long side (depth) of column [m]
L = 175*10^-3; % column length in middle [m]
Lb = 168*10^-3; % column length at bottom [m]
Lt = 75*10^-3; % column length on top [m]

% story heights [m] (from ground to mid floor)
H(1) = Lb + t/2;
for i = 2:5
    H(i) = H(i-1) + L + t;
end

x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:length(omegas)
    subplot(1,length(omegas),i)
    hold on
    plot(phi(:,i),x,'-m')
    if phi(2,i)*SSIphi(1,i) < 0 % Swap sign on mode shape
        plot([0  ;-SSIphi(:,i)],x,'go-.');
        plot(-SSIphi(1:end,i),x(2:end),'g.','markersize',30)
    else
        plot([0  ;SSIphi(:,i)],x,'go-.');
        plot(SSIphi(1:end,i),x(2:end),'g.','markersize',30)
    end
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1 
        legend('Numerical','SSI','Location','northwest')
    end
end

sgtitle('Numerical mode shapes, calibrated by SSI','FontSize',20)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

save('.\datam\Eigenvalue_residual.mat','Knew');
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
xlim([10^-0 10^2+10^2.3])
ylim([10^-0 10^4])
xlabel('norm (Residual)','FontSize',14)
ylabel('norm (Stiffness Change)','FontSize',14)
title('L-curve','FontSize',20)
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