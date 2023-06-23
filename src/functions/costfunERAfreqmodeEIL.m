function J=costfunERAfreqmodeEIL(EIL)
EI = EIL(1);
L = EIL(2);
% Natural frequencies and normalized mode shapes from relevant OMA method 
% ERA
data = load('ERAmodal.mat');
ERAFreq = data.ERAFreq;
omegaOMA = ERAFreq * 2 * pi;
phiOMA = data.phi_ERA;

% story heights [m] (from ground to mid floor)
h = 1*10^-3; % short side of column [m]
b = 30*10^-3; % long side (depth) of column [m]
t = 0.015; % floor height [m]
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

%%%%%%%%%%%%%%% hvad sker der her %%%%%%%%%%%%%%%%%%%
% gravitational force on each floor [N] 
P = flip(cumsum(m))*g;

k0 = sqrt(P./(EI)); % parameter k in DE [1/m]
F = 1; % imposed horizontal load [N] 
% constants boundary conditions for a cantilever beam (homogen solution)
c4 = F./(EI.*k0.^3); % randbetingelse for forskydning w'''(L)=F
c3 = c4.*(cos(k0.*Lh)-1)./sin(k0.*Lh); % randbetingelse for ingen moment w''(L)=0 
%c3 = -c4 * sin(k0.*L)/cos(k0.*L); % randbetingelse for ingen moment w''(L)=0 
c2 = -c4; % randbetingelse for ingen rotation w'(0)=0 
c1 = -c3; % randbetingelse for ingen flytning w(0)=0 
% deflection from imposed load
wL = c1 + c2.*k0.*Lh + c3.*cos(k0.*Lh) + c4.*sin(k0.*Lh); %[m]
% lateral stiffness
k = F./wL; % [N/m]

% stiffness matrix
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);
% mass matrix 
M = m.*eye(5);
% Mass matrix
% M = [2.3553        0         0         0         0
%          0    2.3690         0         0         0
%          0         0    2.3690         0         0
%          0         0         0    2.3690         0
%          0         0         0         0    2.4467];
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
J=sum((omegas-omegaOMA).^2)*1+sum(sum((abs(U)-abs(phiOMA)).^2))*1;
end