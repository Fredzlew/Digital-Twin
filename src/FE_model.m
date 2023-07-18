% parameters
clc; clear; close all;
addpath(genpath('.\data'))
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

% mass matrix
M = m.*eye(5);

% Constitutive stiffness
kc = 12*EI./Lh.^3;
% Geometric stiffness
kg = 6/5*m*g./Lh;

for i = 1:5
    k(i) = kc(i) - (sum(kg(i+1:5)));
end

% stiffness matrix
for i = 1:4
    K(i,i) = k(i)+k(i+1);
    K(i,i+1) = -k(i+1);
    K(i+1,i) = -k(i+1);
end
K(5,5) = k(5);
disp(K)

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
disp(fn)

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
save('.\data\modelprop.mat','K','k','M','H','U','fn','L','EI');

% plotting the mode shapes
x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:length(omegas)
    subplot(1,length(omegas),i)
    hold on
    plot(phi(:,i),x,'-m')
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
    
end
sgtitle('Numerical (JAN)','FontSize',20) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);

T = array2table([num2cell(x'),num2cell(phi)]);
T.Properties.VariableNames(1:6) = {'height','OMAphi1','OMAphi2','OMAphi3','OMAphi4','OMAphi5'};
writetable(T,'C:\Users\User\Danmarks Tekniske Universitet\Frederik Emil Serritzlew - Kandidat\Data\Kap2_FE_modeshapes.csv')

