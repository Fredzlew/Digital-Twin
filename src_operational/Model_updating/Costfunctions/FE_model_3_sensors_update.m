% parameters
clc; clear; close all;
addpath(genpath('..\..\data'),genpath('.\functions'),genpath('..\..\npy-matlab-master'))
% Loading modal parameters from OMA 
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
x = input(promptt);
setGlobalx(x)
if x == 1
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSIfreq_5_2_1.npy');
    SSIomega = SSIFreq * 2 * pi;
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSImodes_5_2_1.npy');
elseif x == 2
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSIfreq_no_damp.npy');
    SSIomega = SSIFreq * 2 * pi;
    SSIphi = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSImodes_no_damp.npy');
end

for loop = 1:5
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

% local stiffnesses
kc = 12*EI./Lh.^3;
kg = 6/5*m*g./Lh;
for i = 1:5
    k2(i) = kc(i) - (sum(kg(i+1:5)));
end

% global stiffness matrix
for i = 1:4
    Km(i,i) = k2(i)+k2(i+1);
    Km(i,i+1) = -k2(i+1);
    Km(i+1,i) = -k2(i+1);
end
Km(5,5) = k2(5);


if loop == 1
    Stiff = fminsearch(@costfunSSIfreq_3_sensors,k2); % SSI, frequency
elseif loop == 2
    Stiff = fminsearch(@costfunSSImode_3_sensors,k2); % SSI, mode shape
elseif loop == 3
    Stiff = fminsearch(@costfunSSIfreqmode_3_sensors,k2); % SSI, frequency + mode shape
elseif loop == 4
    Stiff = fminsearch(@costfunSSIfreqmodeEIL_3_sensors,EIL); % SSI (JAN), EI + L
elseif loop == 5
    Stiff = fminsearch(@costfunSSImode_mac_3_sensors,k2); % SSI, mode shape
end

if loop == 4
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
    
    % Constitutive stiffness
    kc = 12*EI./Lh.^3;
    % Geometric stiffness
    kg = 6/5*m*g./Lh;
    
    for i = 1:5
        k(i) = kc(i) - (sum(kg(i+1:5)));
    end
elseif (loop == 1) || (loop == 2) || (loop == 3) || (loop == 5) 
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
        Uss(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization
for j = 1:5
    for l = 1:5
        U(l,j) = Uss(l,j)/Uss(5,j);
    end
end % end normalization

% saving the parameters
OMAfreq=SSIFreq;
OMAphi=SSIphi;
if x == 1
if loop == 1
    save('.\data_updated_par_3_sensors\SSIfreq_3_sensors_high.mat','OMAphi','OMAfreq','K','Km','U','fn','k');
elseif loop == 2
    save('.\data_updated_par_3_sensors\SSImode_3_sensors_high.mat','OMAphi','OMAfreq','K','Km','U','fn','k')
elseif loop == 3
    save('.\data_updated_par_3_sensors\SSIfreqmode_3_sensors_high.mat','OMAphi','OMAfreq','K','Km','U','fn','k');
elseif loop == 4
    save('.\data_updated_par_3_sensors\SSIfreqmodeEIL_3_sensors_high.mat','OMAphi','OMAfreq','K','L','EI','Km','U','fn','k');
elseif loop == 5
    save('.\data_updated_par_3_sensors\SSImode_mac_3_sensors_high.mat','OMAphi','OMAfreq','K','Km','U','fn','k')
end
elseif x == 2
if loop == 1
    save('.\data_updated_par_3_sensors\SSIfreq_3_sensors_no_damp.mat','OMAphi','OMAfreq','K','Km','U','fn','k');
elseif loop == 2
    save('.\data_updated_par_3_sensors\SSImode_3_sensors_no_damp.mat','OMAphi','OMAfreq','K','Km','U','fn','k')
elseif loop == 3
    save('.\data_updated_par_3_sensors\SSIfreqmode_3_sensors_no_damp.mat','OMAphi','OMAfreq','K','Km','U','fn','k');
elseif loop == 4
    save('.\data_updated_par_3_sensors\SSIfreqmodeEIL_3_sensors_no_damp.mat','OMAphi','OMAfreq','K','L','EI','Km','U','fn','k');
elseif loop == 5
    save('.\data_updated_par_3_sensors\SSImode_mac_3_sensors_no_damp.mat','OMAphi','OMAfreq','K','Km','U','fn','k')
end
end
end
