%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame eigenvalue and mode shape residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
m = diag(M);
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

% Making the sensitivity matrix

% The sensitivity matrix out from the analytical mode shapes for the eigenvalue:
syms k1 k2 k3 k4 k5 m1 m2 m3 m4 m5

Msym = [m1,0,0,0,0;0,m2,0,0,0;0,0,m3,0,0;0,0,0,m4,0;0,0,0,0,m5];

Ksym = [k1+k2,-k2,0,0,0;-k2,k2+k3,-k3,0,0;0,-k3,k3+k4,-k4,0;0,0,-k4,k4+k5,-k5;0,0,0,-k5,k5];

Gs = zeros(length(omegas),length(k));

% number of modes
nm = length(omegas);
% number of parameters
np = length(k);

pars = [k1,k2,k3,k4,k5];
parm = [m1,m2,m3,m4,m5];
for j = 1:nm
    for kk = 1:np
        par = pars(kk);
        Gs(j,kk) = (U(:,j)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j));
    end
end

Gm = zeros(length(omegas),length(k));
for j = 1:nm
    for kk = 1:np
        par = parm(kk);
        Gm(j,kk) = (U(:,j)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j));
    end
end

Geig = [Gs,Gm];

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

Gs = zeros(25,5);
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
            Gs(ind,kk) = b(1);
            Gs(ind+1,kk) = b(2);
            Gs(ind+2,kk) = b(3);
            Gs(ind+3,kk) = b(4);
            Gs(ind+4,kk) = b(5);
          end
        end
    end
    ind = ind+5;
end



% making the a factor for the mass parameters
for j = 1:nm % looping over the mode shape
    for kk = 1:np % looping over the parameters
        for h = 1:nm % Looping over modes shape included in approximated
             if j == h 
                par = parm(kk);
                am(j,kk,h) = -1/2 * U(:,j)'*(diff(Msym,par))*U(:,j);
            elseif j ~= h 
                par = parm(kk);
                am(j,kk,h) = (U(:,h)'*(-omega_squared(j)*diff(Msym,par)+diff(Ksym,par))*U(:,j))/(omega_squared(j)-omega_squared(h));
            end
           
        end
    end
end

Gm = zeros(25,5);
% now finding the sensivity matrix for the mass parameters
ind = 1;
b1 = zeros(5,1);
b = zeros(5,1);
for j = 1:nm
    for kk = 1:np
        for h = 1:nm
          b1 = am(j,kk,h)*U(:,h);
          b = b+b1;
          if h == nm
            Gm(ind,kk) = b(1);
            Gm(ind+1,kk) = b(2);
            Gm(ind+2,kk) = b(3);
            Gm(ind+3,kk) = b(4);
            Gm(ind+4,kk) = b(5);
          end
        end
    end
    ind = ind+5;
end

% The combined sensivity matrix of mode shape
Gmod = [Gs,Gm];

% The sensivity matrix
G = [Geig;Gmod];

% condition number of G
con = cond(G);

% symmetric weighting matrix
V = ones(size(G,1),1);
%Weps = diag(omegasq_m)^-2;
V(1) = 1;
V(6:end) = 1000;
Weps = eye(size(G,1)).*V;

% a simple version of the parameter weighting matrix
%Wtheta = eye(size(G,2));
% The parameter weighting matrix for regularization problem from equation
% (19) and (20)
Gamma = diag(diag(G'*Weps*G));
Wtheta = mean(diag(Gamma))/mean(diag(Gamma^-1))*Gamma^-1;

% Difne lambda value:
lambda = sqrt(0.3437); 
%lambda = sqrt(12.0331);

% initial parameters changes
dx = zeros(size(G,2),1);

% reshape of modeshape to 25x1
SSIphi_re = reshape(SSIphi,size(SSIphi,1)*size(SSIphi,2),1);

% Iterazation over the parameters
for ii = 1:100
    % updated parameters
    k = k+dx(1:5);
    m = m+dx(6:end);
    
    % stiffness matrix
    for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
    end
    K(5,5) = k(5);
    
    % Mass matrix 

    M = [m(1), 0, 0, 0, 0;...
         0, m(2), 0, 0, 0;...
         0, 0, m(3), 0, 0;...
         0, 0, 0, m(4), 0;...
         0, 0, 0, 0, m(5)];

    % eigenvalue problem
    [U,D] = eig(K,M);
    
    % natural frequencies from eigenvalues
    omega = real(sqrt(diag(D)));
    
    % sort frequencies and mode shapes
    [~,iw] = sort(omega);

    % natural frequencies [rad/s]
    omegas = omega(iw);
    
    % Frequencie
    fn = omegas/(2*pi);

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
        
    % normalization of the residual eig
%     MVec_x = max(reig(:,ii)); % start normalization
%     mVec_x = min(reig(:,ii));
% 
%     if abs(MVec_x(1)) > abs(mVec_x(1))
%         mxVec_x(1) = MVec_x(1);
%     else
%         mxVec_x(1) = mVec_x(1);
%     end
%     for l = 1:length(reig(:,ii))
%         reig(l,ii) = reig(l,ii)/mxVec_x(1);
%     end % end normalization

    % Combines residuals
    r(:,ii) = [reig(:,ii);rmod(:,ii)];

    disp(sum(abs(r(:,ii))))

    
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


% Display accuracy of mode shapes
% CrossMAC plot of mode shapes
for i = 1:100
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
knew = k+dx(1:5);
mnew = m+dx(6:end);

% Stiffness matrix
for i = 1:4
    Knew(i,i) = knew(i)+knew(i+1);
    Knew(i,i+1) = -knew(i+1);
    Knew(i+1,i) = -knew(i+1);
end
Knew(5,5) = k(5);

Mnew = [mnew(1), 0, 0, 0, 0;...
        0, mnew(2), 0, 0, 0;...
        0, 0, mnew(3), 0, 0;...
        0, 0, 0, mnew(4), 0;...
        0, 0, 0, 0, mnew(5)];

% eigenvalue problem
[U,D] = eig(Knew,Mnew);

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

sgtitle('Mode shapes','FontSize',20)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% Display accuracy of frequency
disp(strcat('Frequency accuracy,1 : ',num2str(min(SSIFreq(1),fn(1))/max(SSIFreq(1),fn(1))*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(min(SSIFreq(2),fn(2))/max(SSIFreq(2),fn(2))*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(min(SSIFreq(3),fn(3))/max(SSIFreq(3),fn(3))*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(min(SSIFreq(4),fn(4))/max(SSIFreq(4),fn(4))*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(min(SSIFreq(5),fn(5))/max(SSIFreq(5),fn(5))*100),'%'));
disp(strcat('Mean frequency accuracy : ',num2str(mean(min(SSIFreq,fn)./max(SSIFreq,fn)*100)),'%'));

% Display accuracy of mode shapes
% CrossMAC plot of mode shapes
mac=crossMAC(U,SSIphi,1,[SSIFreq,fn]);
dmac = diag(mac);
disp('Modal Assurance Criterion between Numerical modeshapes and SSI  : ')
disp(strcat(num2str(mac)));
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy (MAC),1 : ',num2str(dmac(1)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),2 : ',num2str(dmac(2)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),3 : ',num2str(dmac(3)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),4 : ',num2str(dmac(4)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),5 : ',num2str(dmac(5)*100),'%'));
disp(strcat('Mean mode shape accuracy (MAC): ',num2str(mean(dmac)*100),'%'));
disp('----------------------------------------------------------------------')

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
% xlim([10^-4 10^-1+10^-1/2])
% ylim([10^-3 10^-0+10^-0/1.5])
xlabel('norm (Residual)')
ylabel('norm (Stiffness Change)')
title('L-curve')
% Finding the optimal value for lambda
Val = 102.716;
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