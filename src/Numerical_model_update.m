% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
%    column with lumbed masses     
%         
%        o <-- H(5),m(5)
%        |
%        o <-- H(4),m(4)
%   y    |
%   Ʌ    o <-- H(3),m(3)
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
FDDphi = FDD.phi_FDD;


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
EIL = [EI,L];

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
% Initial stiffness matrix
% stiffness matrix
for i = 1:4
    Km(i,i) = k2(i)+k2(i+1);
    Km(i,i+1) = -k2(i+1);
    Km(i+1,i) = -k2(i+1);
end
Km(5,5) = k2(5);

% Change cost function to use correct natural frequencies
% Define what OMA method is used (also change data in costfunction)
% MODE = 2; % 1=SSI, 2=ERA, 3=FDD
prompt = "Which OMA is used (1=SSI, 2=ERA, 3=FDD)? ";
MODE = input(prompt);
% Define and minimize cost function
promptt = "Which costfun is used to find calibrated stiffness (1=SSIFreq, 2=ERAFreq, 3=FDDFreq, 4=SSImodes, 5=ERAmodes, 6=FDDmodes, 7=SSIFreqmodes, 8=ERAFreqmodes, 9=FDDFreqmodes )? ";
Stif = input(promptt);
if Stif == 1
    Stiff = fminsearch(@costfunSSIfreq,k2); % SSI, frequency
elseif Stif == 2
    Stiff = fminsearch(@costfunERAfreq,k2); % ERA, frequency
elseif Stif == 3
    Stiff = fminsearch(@costfunFDDfreq,k2); % FDD, frequency
elseif Stif == 4
    Stiff = fminsearch(@costfunSSImode,k2); % SSI, mode shape
elseif Stif == 5
    Stiff = fminsearch(@costfunERAmode,k2); % ERA, mode shape
elseif Stif == 6
    Stiff = fminsearch(@costfunFDDmode,k2); % FDD, mode shape
elseif Stif == 7
    Stiff = fminsearch(@costfunSSIfreqmode,k2); % SSI, frequency + mode shape
elseif Stif == 8
    Stiff = fminsearch(@costfunERAfreqmode,k2); % ERA, frequency + mode shape
elseif Stif == 9
    Stiff = fminsearch(@costfunFDDfreqmode,k2); % FDD, frequency + mode shape
elseif Stif == 10
    Stiff = fminsearch(@costfunSSIfreqmodeEIL,EIL); % SSI, EI + L
elseif Stif == 11
    Stiff = fminsearch(@costfunERAfreqmodeEIL,EIL); % ERA, EI + L
elseif Stif == 12
    Stiff = fminsearch(@costfunFDDfreqmodeEIL,EIL); % FDD, EI + L
elseif Stif == 13
    Stiff = fminsearch(@costfunSSIfreqmodeEILJAN,EIL); % SSI (JAN), EI + L
elseif Stif == 14
    Stiff = fminsearch(@costfunERAfreqmodeEILJAN,EIL); % ERA (JAN), EI + L
elseif Stif == 15
    Stiff = fminsearch(@costfunFDDfreqmodeEILJAN,EIL); % FDD (JAN), EI + L
end

if (Stif == 10) || (Stif == 11) || (Stif == 12)
    EI = Stiff(1);
    L = Stiff(2);
    Lb = 168/175*L; % column length at bottom [m]
    Lt = 75/175*L; % column length on top [m]

    % story heights [m] (from ground to mid floor)
    H(1) = Lb + t/2;
    for i = 2:5
        H(i) = H(i-1) + L + t;
    end

    % element properties
    g = 9.81; % [m/s^2]
    
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
    k = F./wL; % [N/m]
elseif (Stif == 13) || (Stif == 14) || (Stif == 15)
    EI = Stiff(1);
    L = Stiff(2);
    Lb = 168/175*L; % column length at bottom [m]
    Lt = 75/175*L; % column length on top [m]

    % story heights [m] (from ground to mid floor)
    H(1) = Lb + t/2;
    for i = 2:5
        H(i) = H(i-1) + L + t;
    end

    % element properties
    g = 9.81; % [m/s^2]
    
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
    m = mf+4*b*h*rho*[Lh(1) Lh(2) Lh(3) Lh(4) Lh(5)+Lt+t/2]; 
    
    kc = 12*EI/L^3;
    kg = 6/5*m*g/L;
    for i = 1:5
        k(i) = kc - (sum(kg(i+1:5)));
    end
elseif (Stif == 1) || (Stif == 2) || (Stif == 3) || (Stif == 4) || (Stif == 5) || (Stif == 6) || (Stif == 7) || (Stif == 8) || (Stif == 9) 
    % Define optimal stiffness
    k = Stiff;
end

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
    if phi(2,i)*OMAphi(1,i) < 0 % Swap sign on mode shape
        plot([0  ;-OMAphi(:,i)],x,'go-.');
        plot(-OMAphi(1:end,i),x(2:end),'g.','markersize',30)
    else
        plot([0  ;OMAphi(:,i)],x,'go-.');
        plot(OMAphi(1:end,i),x(2:end),'g.','markersize',30)
    end
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
    OMAphi=SSIphi;
elseif MODE==2
    OMAfreq=ERAFreq;
    OMAphi=ERAphi;
else
    OMAfreq=FDDFreq;
    OMAphi=FDDphi;
end
disp(strcat('Frequency accuracy,1 : ',num2str(min(OMAfreq(1),fn(1))/max(OMAfreq(1),fn(1))*100),'%'));
disp(strcat('Frequency accuracy,2 : ',num2str(min(OMAfreq(2),fn(2))/max(OMAfreq(2),fn(2))*100),'%'));
disp(strcat('Frequency accuracy,3 : ',num2str(min(OMAfreq(3),fn(3))/max(OMAfreq(3),fn(3))*100),'%'));
disp(strcat('Frequency accuracy,4 : ',num2str(min(OMAfreq(4),fn(4))/max(OMAfreq(4),fn(4))*100),'%'));
disp(strcat('Frequency accuracy,5 : ',num2str(min(OMAfreq(5),fn(5))/max(OMAfreq(5),fn(5))*100),'%'));
disp(strcat('Mean frequency accuracy : ',num2str(mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)),'%'));

% Display accuracy of mode shapes
[MSacc,TOTacc]=modeshapeacc(OMAphi,U);
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy,1 : ',num2str(MSacc(1)*100),'%'));
disp(strcat('Mode shape accuracy,2 : ',num2str(MSacc(2)*100),'%'));
disp(strcat('Mode shape accuracy,3 : ',num2str(MSacc(3)*100),'%'));
disp(strcat('Mode shape accuracy,4 : ',num2str(MSacc(4)*100),'%'));
disp(strcat('Mode shape accuracy,5 : ',num2str(MSacc(5)*100),'%'));
disp(strcat('Mean mode shape accuracy : ',num2str(TOTacc*100),'%'));


% Display stiffness matrices and changes in these
disp('----------------------------------------------------------------------')
disp('Original stiffness matrix [N/m]  : ')
disp(strcat(num2str(Km)));
disp('Calibrated stiffness matrix [N/m]  : ')
disp(strcat(num2str(K)));
disp('Changes in stiffness matrix [N/m] : ')
disp(strcat(num2str(abs(abs(K)-abs(Km)))));
disp('Total change in stiffness matrix [N/m] : ')
disp(strcat(num2str(sum(sum(abs(abs(K)-abs(Km)))))));
disp('----------------------------------------------------------------------')
disp(strcat('Total mean accuracy (mode+freq) : ',num2str(mean([TOTacc*100,mean(min(OMAfreq,fn)./max(OMAfreq,fn)*100)])),'%'));
disp('----------------------------------------------------------------------')

% Histogram of the different stiffness from different OMA methods and cost
% functions
% Stiffnesses
K_SSI_freq = [6259.3525,6182.6736,6679.9769,7579.0717,4016.2825];
K_ERA_freq = [6710.4572,5363.6827,5569.9894,8039.6749,4816.2026];
K_FDD_freq = [6322.7709,5997.1983,6490.1024,7738.5212,4168.6269];
K_SSI_mode = [6583.7572,6437.9312,6684.1129,7230.6171,3830.8476];
K_ERA_mode = [6557.7385,6419.7800,6640.9573,7223.8272,3847.3987];
K_FDD_mode = [6677.7406,6432.6900,6642.5802,7227.1624,3869.3218];
K_SSI_freq_mode = [6344.3278,6431.9967,6852.2047,7311.6168,3768.8334];
K_ERA_freq_mode = [6189.0799,6341.3222,6853.1140,7358.7121,3725.9845];
K_FDD_freq_mode = [6347.0535,6409.5619,6836.3663,7333.9437,3777.6936];
y = [K_SSI_freq;K_ERA_freq;K_FDD_freq;K_SSI_mode;K_ERA_mode;K_FDD_mode;K_SSI_freq_mode;K_ERA_freq_mode;K_FDD_freq_mode;k2];
% Find value at each tip
figure
hold on
b = bar(y,'stacked');
for i = 1:size(y,2)
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = string(b(i).YData);
    text(xtips,ytips./1.1,labels,'HorizontalAlignment','center','VerticalAlignment','middle')
end
hold off
legend('k1','k2','k3','k4','k5')
title('Values of stiffness for different OMA methods and cost functions')
xlabel('Method')
xticks(xtips)
xticklabels({'SSI (freq)','ERA (freq)','FDD (freq)','SSI (mode)','ERA (mode)','FDD (mode)','SSI (freq+mode)','ERA (freq+mode)','FDD (freq+mode)','Geometric Stiffness'})
ylabel('Stiffness [N/m]')

% 3D plot
figure
bar3(y)
title('Values of stiffness for different OMA methods and cost functions')
%xlabel('ks')
%ylabel('Method')
xticks([1,2,3,4,5])
xticklabels({'k1','k2','k3','k4','k5'})
yticks(xtips)
yticklabels({'SSI (freq)','ERA (freq)','FDD (freq)','SSI (mode)','ERA (mode)','FDD (mode)','SSI (freq+mode)','ERA (freq+mode)','FDD (freq+mode)','Geometric Stiffness'})
zlabel('Stiffness [N/m]')

% CrossMAC plot of mode shapes
mac=crossMAC(U,OMAphi,MODE,[OMAfreq,fn]);
dmac = diag(mac);
if MODE==1
    disp('Modal Assurance Criterion between Numerical modeshapes and SSI  : ')
    disp(strcat(num2str(mac)));
elseif MODE==2
    disp('Modal Assurance Criterion between Numerical modeshapes and ERA  : ')
    disp(strcat(num2str(mac)));
else
    disp('Modal Assurance Criterion between Numerical modeshapes and FDD  :' )
    disp(strcat(num2str(mac)));
end
disp('----------------------------------------------------------------------')
disp(strcat('Mode shape accuracy (MAC),1 : ',num2str(dmac(1)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),2 : ',num2str(dmac(2)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),3 : ',num2str(dmac(3)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),4 : ',num2str(dmac(4)*100),'%'));
disp(strcat('Mode shape accuracy (MAC),5 : ',num2str(dmac(5)*100),'%'));
disp(strcat('Mean mode shape accuracy (MAC): ',num2str(mean(dmac)*100),'%'));
disp('----------------------------------------------------------------------')

if Stif == 1
    save('.\data\costfunupdateSSIfreq.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 2
    save('.\data\costfunupdateERAfreq.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 3
    save('.\data\costfunupdateFDDfreq.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 4
    save('.\data\costfunupdateSSImode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 5
    save('.\data\costfunupdateERAmode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 6
    save('.\data\costfunupdateFDDmode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 7
    save('.\data\costfunupdateSSIfreqmode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 8
    save('.\data\costfunupdateERAfreqmode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 9
    save('.\data\costfunupdateFDDfreqmode.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 10
    save('.\data\costfunupdateSSIfreqmodeEIL.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 11
    save('.\data\costfunupdateERAfreqmodeEIL.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 12
    save('.\data\costfunupdateFDDfreqmodeEIL.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 13
    save('.\data\costfunupdateSSIfreqmodeEILJAN.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 14
    save('.\data\costfunupdateERAfreqmodeEILJAN.mat','OMAphi','OMAfreq','K','x','phi');
elseif Stif == 15
    save('.\data\costfunupdateFDDfreqmodeEILJAN.mat','OMAphi','OMAfreq','K','x','phi');
end
