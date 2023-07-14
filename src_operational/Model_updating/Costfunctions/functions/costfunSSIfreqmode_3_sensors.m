function J=costfunSSIfreqmode(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
f = getGlobalx;
% Natural frequencies and normalized mode shapes from relevant OMA method 
% SSI
if f == 1
    % High damping
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSIfreq_5_2_1.npy');
    omegaOMA = SSIFreq * 2 * pi;
    phiOMA = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSImodes_5_2_1.npy');
elseif f == 2
    % No damping
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSIfreq_no_damp.npy');
    omegaOMA = SSIFreq * 2 * pi;
    phiOMA = readNPY('..\..\data\experimental_data\Modal_par_3_sensors\SSImodes_no_damp.npy');
end
% Stiffness matrix
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
        Uss(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization
for j = 1:5
    for l = 1:5
        U(l,j) = Uss(l,j)/Uss(5,j);
    end
end % end normalization
U_3_sensors=[U(2,:);U(3,:);U(5,:)];
J=sum((omegas-omegaOMA).^2)*1+sum(sum((abs(U_3_sensors)-abs(phiOMA)).^2))*1;
end