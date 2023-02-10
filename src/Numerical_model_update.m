% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%    column with lumbed masses     
%         
%        o <-- H(5),m(5)
%        |
%        o <-- H(4),m(4)
%        |
%   y    o <-- H(3),m(3)
%  / \   |
%   |    o <-- H(2),m(2)
%   |    |
%   |    o <-- H(1),m(1)
%   |   _|_    ------> x
%      //// 
% ----------------------------------------------------
% Loading modal parameters from OMA 
SSI = load('SSImodal.mat'); % loading modal parameters from SSI
ERA = load('ERAmodal.mat'); % loading modal parameters from ERA
FDD = load('FDDmodal.mat'); % loading modal parameters from FDD

SSIFreq = SSI.SSIFreq;
SSIomega = SSIFreq * 2 * pi;
SSIphi = SSI.phi_SSI;

ERAFreq = ERA.ERAFreq;
ERAomega = ERAFreq * 2 * pi;
ERAphi = ERA.phi_ERA;

FDDFreq = FDD.fn';
FDDomega = FDDFreq * 2 * pi;
FDDphi = FDD.phi_FDD';


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

% element properties
g = 9.81; % [m/s^2]
I = 4*1/12*b*h^3; % moment of inertia of 4 columns [m^4]
E = 210*10^9; % [Pa]
EI = E*I; % [N*m^2]

% storage height for each floor [m]
L(1) = H(1);
for i = 2:5
    L(i) = H(i)-H(i-1);
end

% lumbed masses [kg]
ml = 0.134; % list
mp = 1.906; % plate
mb = 0.0160; % bolts
mf = mp + 2 * ml + mb; % total mass of 1 floor (plate + 2 lists and bolts)
% mf = 2.19; % floor (plate + 2 lists and bolts)
rho = 7850; % density column [kg/m^3]

% total mass of frame [kg]
m = mf+4*b*h*rho*[L(1) L(2) L(3) L(4) L(5)+Lt+t/2]; 

%%%%%%%%%%%%%%% hvad sker der her %%%%%%%%%%%%%%%%%%%
% gravitational force on each floor [N] 
P = flip(cumsum(m))*g;

k0 = sqrt(P./(EI)); % parameter k in DE [1/m]
F = 1; % imposed horizontal load [N] 
% constants boundary conditions for a cantilever beam (homogen solution)
c4 = F./(EI.*k0.^3); % randbetingelse for forskydning w'''(L)=F
c3 = c4.*(cos(k0.*L)-1)./sin(k0.*L); % randbetingelse for ingen moment w''(L)=0 
%c3 = -c4 * sin(k0.*L)/cos(k0.*L); % randbetingelse for ingen moment w''(L)=0 
c2 = -c4; % randbetingelse for ingen rotation w'(0)=0 
c1 = -c3; % randbetingelse for ingen flytning w(0)=0 
% deflection from imposed load
wL = c1 + c2.*k0.*L + c3.*cos(k0.*L) + c4.*sin(k0.*L); %[m]
% lateral stiffness
k2 = F./wL; % [N/m]

% Change cost function to use correct natural frequencies
% Define what OMA method is used (also change data in costfunction)
MODE = 1; % 1=SSI, 2=ERA, 3=FDD

% Define and minimize cost function
Stiff = fminsearch(@costfunOMA,k2);
% Define optimal stiffnesses
k = Stiff;
% stiffness matrix
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);

% mass matrix 
M = m.*eye(5);

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
% frequencies [Hz]
fn = omegas./(2*pi);

% normalization
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        U(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization

%... Save results in *.mat file .................
if MODE==1
    save('.\data\modelpropupdateSSI.mat','K','M','H','U');
elseif MODE==2
    save('.\data\modelpropupdateERA.mat','K','M','H','U');
else
    save('.\data\modelpropupdateFDD.mat','K','M','H','U');
end

% plotting the mode shapes
x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];
if MODE==1
    OMAphi = SSIphi;
elseif MODE==2
    OMAphi = ERAphi;
else
    OMAphi = FDDphi;
end
for i=1:length(omegas)
    subplot(1,length(omegas),i)
    hold on
    plot(phi(:,i),x,'-m')
    plot([0  ;OMAphi(:,i)],x,'go-.');
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    if i==1 && MODE==1
        legend('Numerical','SSI','Location','northwest')
    elseif i==1 && MODE==2
        legend('Numerical','ERA','Location','northwest')
    elseif i==1
        legend('Numerical','FDD','Location','northwest')
    end
end
if MODE==1
    sgtitle('Numerical mode shapes, calibrated by SSI','FontSize',20)
elseif MODE==2
    sgtitle('Numerical mode shapes, calibrated by ERA','FontSize',20)
else
    sgtitle('Numerical mode shapes, calibrated by FDD','FontSize',20)
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

% Display accuracy of natural frequencies
if MODE==1
    OMAfreq=SSIFreq;
elseif MODE==2
    OMAfreq=ERAFreq;
else
    OMAfreq=FDDFreq;
end
disp(strcat('Frequency accuracy,1 : ',num2str(OMAfreq(1)/fn(1)*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(OMAfreq(2)/fn(2)*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(OMAfreq(3)/fn(3)*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(OMAfreq(4)/fn(4)*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(OMAfreq(5)/fn(5)*100),'%'));

% Display accuracy of mode shapes
% disp(strcat('Mode shape accuracy,1 : ',num2str(mean(abs(SSIphi(:,1))./abs(U(:,1)))*100),'%'));
% disp(strcat('Mode shape accuracy,2 : ',num2str(mean(abs(SSIphi(:,2))./abs(U(:,2)))*100),'%'));
% disp(strcat('Mode shape accuracy,3 : ',num2str(mean(abs(SSIphi(:,3))./abs(U(:,3)))*100),'%'));
% disp(strcat('Mode shape accuracy,4 : ',num2str(mean(abs(SSIphi(:,4))./abs(U(:,4)))*100),'%'));
% disp(strcat('Mode shape accuracy,5 : ',num2str(mean(abs(SSIphi(:,5))./abs(U(:,5)))*100),'%'));


