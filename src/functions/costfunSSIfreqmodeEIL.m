function J=costfunSSIfreqmodeEIL(EI,L)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
% Normal force on bottom floor [total weight in kg]
N = -2.4467*9.82;     % minus = compression
% Natural frequencies and normalized mode shapes from relevant OMA method 
% SSI
omegaOMA =  [10.3880;
   31.5746;
   49.6714;
   63.5516;
   72.9177];
phiOMA =  [0.2939    0.7129    1.0000   -0.9947    0.4457
    0.6249    1.0000    0.2185    1.0000   -0.8244
    0.7927    0.4622   -0.8788    0.0171    1.0000
    0.9963   -0.3065   -0.4272   -0.9361   -0.9710
    1.0000   -0.7890    0.6735    0.5416    0.3842];
% Constitutive stiffness matrix
Ke = [E*A/L     0           0           -E*A/L  0           0
          0     12*E*I/L^3  6*E*I/L^2   0       -12*E*I/L^3 6*E*I/L^2
          0     6*E*I/L^2   4*E*I/L     0       -6*E*I/L^2  2*E*I/L
     -E*A/L     0           0           E*A/L   0           0
          0     -12*E*I/L^3 -6*E*I/L^2  0       12*E*I/L^3  -6*E*I/L^2
          0     6*E*I/L^2   2*E*I/L     0       -6*E*I/L^2  4*E*I/L     ];
% Geometrical stiffness matrix 
Kg =(N/(30*L))*[ 36   3*L   -36   3*L
      3*L  4*L^2 -3*L -L^2   
     -36  -3*L   36   -3*L   
      3*L -L^2  -3*L  4*L^2 ];
% Total stiffness matrix
K = Ke+Kg;
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