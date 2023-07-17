%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5 storey frame eigenvalue and mode shape residual %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath(genpath('..\..\npy-matlab-master'),genpath('..\..\data'),genpath('..\Costfunctions\functions'))
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

% Height of building
H = filename.H;

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
    
    % Frequency [Hz]
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


    % plotting the mode shapes
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
        title(['Freq. Approx = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        subtitle(['(FE freq. = ' num2str(SSIFreq(i)) ' Hz)'])
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])
        if i==1 
            legend('Numerical','SSI','Location','northwest')
        end
    end
    sgtitle(['Iteration: ' num2str(ii)],'FontSize',20,'FontWeight','Bold');
    if q == 1
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_high_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_high_' num2str(ii) '.png']);
    elseif q == 2
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_no_damp_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_no_damp_' num2str(ii) '.png']);
    end
    close
end
% requirement
req = rank(G'*Weps*G) == size(G,2);
disp(req);

% Convergence plot for the relative failure for frequencies
for ii = 1:ni
    figure (1)
    colororder(gcf,[0.0,0.4470,0.7410;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560;0.4660,0.6740,0.1880])
    hold on
    for j = 1:5
        err(j,ii) = abs(r(j,ii))/omegasq_m(j);
    end
    plot(err(1,1:ii)')%,'r')
    plot(err(2,1:ii)')%,'g')
    plot(err(3,1:ii)')%,'b')
    plot(err(4,1:ii)')%,'c')
    plot(err(5,1:ii)')%,'m')
    xlim([0 100])
    ylim([0 0.1])
    title(['Iteration: ' num2str(ii)],'FontSize',20,'FontWeight','Bold');
    legend('Frequency 1', 'Frequency 2','Frequency 3', 'Frequency 4','Frequency 5')
    xlabel('Iterations [-]')
    ylabel('Relative error [-]')
    hold off
    if q == 1
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_high_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_error_high_' num2str(ii) '.png']);
    elseif q == 2
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_no_damp_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_error_low_' num2str(ii) '.png']);
    end
end


% Display accuracy of mode shapes
% CrossMAC plot of mode shapes
for i = 1:ni
    mac=crossMAC(SSIphi,Umac(:,:,i));
    for j = 1:5
        dmac = diag(mac);
        acc(j,i) = dmac(j);
    end
end

% Convergence plot for the relative failure for mode shapes
for ii = 1:ni
    figure (2)
    colororder(gcf,[0.0,0.4470,0.7410;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560;0.4660,0.6740,0.1880])
    hold on
    plot(acc(1,1:ii)')%,'r')
    plot(acc(2,1:ii)')%,'g')
    plot(acc(3,1:ii)')%,'b')
    plot(acc(4,1:ii)')%,'c')
    plot(acc(5,1:ii)')%,'m')
    xlim([0 100])
    ylim([0.98 1])
    title(['Iteration: ' num2str(ii)],'FontSize',20,'FontWeight','Bold');
    legend('Mode shape 1', 'Mode shape 2','Mode shape 3', 'Mode shape 4','Mode shape 5','Location','southwest')
    xlabel('Iterations [-]')
    ylabel('MAC [-]')
    hold off
    if q == 1
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_high_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_mac_high_' num2str(ii) '.png']);
    elseif q == 2
        %saveas(gcf,['C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\GIF\GIF_no_damp_' num2str(ii) '.png']);
        saveas(gcf,['C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\GIF\GIF_mac_low_' num2str(ii) '.png']);
    end
end

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
    save('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_high.mat','knew','Knew','U','fn','acc','err');
elseif q == 2
    save('.\data_updated_par_sens\Eigenvalue_Mode_shape_residual_no_damp.mat','knew','Knew','U','fn','acc','err');
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
