% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'),genpath('python'))
for loop =1:10
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
SSIFreq = readNPY('omega.npy');
SSIomega = SSIFreq * 2 * pi;
SSIphi = readNPY('phi.npy');

FDDFreq = readNPY('FDDomega.npy');
FDDomega = FDDFreq * 2 * pi;
FDDphi = readNPY('FDDphi.npy');

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
Lh(1) = H(1);
for i = 2:5
    Lh(i) = H(i)-H(i-1);
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

% prompttt = "Using analystic or Total inter-storey stiffness: (1=ANALYSTIC, 2=INTER-STOREY)? ";
% prop = input(prompttt);
% if prop == 1
% %%%%%%%%%%%%%%% hvad sker der her %%%%%%%%%%%%%%%%%%%
% % gravitational force on each floor [N] 
% P = flip(cumsum(m))*g;
% 
% k0 = sqrt(P./(EI)); % parameter k in DE [1/m]
% F = 1; % imposed horizontal load [N] 
% % constants boundary conditions for a cantilever beam (homogen solution)
% c4 = F./(EI.*k0.^3); % randbetingelse for forskydning w'''(L)=F
% c3 = c4.*(cos(k0.*Lh)-1)./sin(k0.*Lh); % randbetingelse for ingen moment w''(L)=0 
% %c3 = -c4 * sin(k0.*L)/cos(k0.*L); % randbetingelse for ingen moment w''(L)=0 
% c2 = -c4; % randbetingelse for ingen rotation w'(0)=0 
% c1 = -c3; % randbetingelse for ingen flytning w(0)=0 
% % deflection from imposed load
% wL = c1 + c2.*k0.*Lh + c3.*cos(k0.*Lh) + c4.*sin(k0.*Lh); %[m]
% % lateral stiffness
% k2 = F./wL; % [N/m]
% % Initial stiffness matrix
% % stiffness matrix
% elseif prop == 2
    kc = 12*EI/L^3;
    kg = 6/5*m*g/L;
    for i = 1:5
        k2(i) = kc - (sum(kg(i+1:5)));
    end
% end
for i = 1:4
    Km(i,i) = k2(i)+k2(i+1);
    Km(i,i+1) = -k2(i+1);
    Km(i+1,i) = -k2(i+1);
end
Km(5,5) = k2(5);
stivhed = diag(Km);
% Change cost function to use correct natural frequencies
% Define what OMA method is used (also change data in costfunction)
% MODE = 2; % 1=SSI, 2=ERA, 3=FDD
% prompt = "Which OMA is used (1=SSI, 3=FDD)? ";
% MODE = input(prompt);
if  (loop == 1) || (loop == 3) || (loop == 5) || (loop == 7) || (loop == 9) 
    MODE = 1;
else
    MODE = 3;
end


% Define and minimize cost function
% promptt = "Which costfun is used to find calibrated stiffness (1=SSIFreq, 2=ERAFreq, 3=FDDFreq, 4=SSImodes, 5=ERAmodes, 6=FDDmodes, 7=SSIFreqmodes, 8=ERAFreqmodes, 9=FDDFreqmodes )? ";
% Stif = input(promptt);
Stif = loop;
if Stif == 1
    Stiff = fminsearch(@costfunSSIfreq,k2); % SSI, frequency
elseif Stif == 2
    Stiff = fminsearch(@costfunFDDfreq,k2); % FDD, frequency
elseif Stif == 3
    Stiff = fminsearch(@costfunSSImode,k2); % SSI, mode shape
elseif Stif == 4
    Stiff = fminsearch(@costfunFDDmode,k2); % FDD, mode shape
elseif Stif == 5
    Stiff = fminsearch(@costfunSSIfreqmode,k2); % SSI, frequency + mode shape
elseif Stif == 6
    Stiff = fminsearch(@costfunFDDfreqmode,k2); % FDD, frequency + mode shape
elseif Stif == 7
    Stiff = fminsearch(@costfunSSIfreqmodeEILJAN,EIL); % SSI (JAN), EI + L
elseif Stif == 8
    Stiff = fminsearch(@costfunFDDfreqmodeEILJAN,EIL); % FDD (JAN), EI + L
end

if (Stif == 7) || (Stif == 8) 
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
    Lh(1) = H(1);
    for i = 2:5
        Lh(i) = H(i)-H(i-1);
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
elseif (Stif == 1) || (Stif == 2) || (Stif == 3) || (Stif == 4) || (Stif == 5) || (Stif == 6) 
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


% plotting the mode shapes
x = [0, H];
phi = [zeros(1,length(U)); U];

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

if Stif == 1
    save('.\data\costfunupdateSSIfreq.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 2
    save('.\data\costfunupdateFDDfreq.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 3
    save('.\data\costfunupdateSSImode.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn')
elseif Stif == 4
    save('.\data\costfunupdateFDDmode.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 5
    save('.\data\costfunupdateSSIfreqmode.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 6
    save('.\data\costfunupdateFDDfreqmode.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 7
    save('.\data\costfunupdateSSIfreqmodeEILJAN.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
elseif Stif == 8
    save('.\data\costfunupdateFDDfreqmodeEILJAN.mat','OMAphi','OMAfreq','K','Km','x','phi','stivhed','H','U','fn');
end
end