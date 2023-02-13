function J=costfunERAfreqmode(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
% Natural frequencies and normalized mode shapes from relevant OMA method 
% ERA
omegaOMA =  [10.0667;
   31.4607;
   49.1633;
   63.2812;
   72.8966];
phiOMA =  [0.2932    0.7087    1.0000   -0.9625    0.4548
    0.6259    1.0000    0.2209    1.0000   -0.7618
    0.7899    0.4613   -0.8978    0.0153    1.0000
    0.9957   -0.3049   -0.4457   -0.9599   -0.9410
    1.0000   -0.7915    0.6830    0.5690    0.4714];
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
        U(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization
J=sum((omegas-omegaOMA).^2)*0.5+sum(sum((abs(U)-abs(phiOMA)).^2))*0.5;
end