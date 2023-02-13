function J=costfunFDDmode(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
% Normalized mode shapes from relevant OMA method 
% FDD
phiOMA =  [0.3034    0.7332    1.0000    1.0000    0.4855
    0.6137    1.0000    0.2127   -0.9640   -0.8445
    0.7760    0.4656   -0.8404   -0.0218    1.0000
    0.9744   -0.3042   -0.4059    0.9037   -0.9874
    1.0000   -0.8046    0.6585   -0.5401    0.3962];
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
J=sum(sum((abs(U)-abs(phiOMA)).^2));
end