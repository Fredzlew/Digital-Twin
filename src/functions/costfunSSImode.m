function J=costfunSSImode(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
% Normalized mode shapes from relevant OMA method 
% SSI
phiOMA =  [0.2939    0.7129    1.0000   -0.9947    0.4457
    0.6249    1.0000    0.2185    1.0000   -0.8244
    0.7927    0.4622   -0.8788    0.0171    1.0000
    0.9963   -0.3065   -0.4272   -0.9361   -0.9710
    1.0000   -0.7890    0.6735    0.5416    0.3842];
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
J=sum(sum((U-phiOMA).^2));
end