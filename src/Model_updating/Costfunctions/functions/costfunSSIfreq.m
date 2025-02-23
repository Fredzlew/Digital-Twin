function J=costfunSSIfreq(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
f = getGlobalx;
% Natural frequencies from relevant OMA method 
% SSI
if f == 1
    % High damping
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_5_2_1.npy');
    omegaOMA = SSIFreq * 2 * pi;
    omegaOMAsq = omegaOMA.^2;
elseif f == 2
    % No damping
    SSIFreq = readNPY('..\..\data\experimental_data\Modal_par\SSIfreq_no_damp.npy');
    omegaOMA = SSIFreq * 2 * pi;
    omegaOMAsq = omegaOMA.^2;
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
omegassq = omegas.^2;
J=sum((omegassq-omegaOMAsq).^2);
end