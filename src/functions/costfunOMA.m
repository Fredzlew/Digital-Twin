function J=costfunOMA(k)
% Mass matrix
M = [2.3553        0         0         0         0
         0    2.3690         0         0         0
         0         0    2.3690         0         0
         0         0         0    2.3690         0
         0         0         0         0    2.4467];
% Natural frequencies from relevant OMA method 
% (SSI)
omegaOMA =  [10.3880;
   31.5746;
   49.6714;
   63.5516;
   72.9177];
% ERA
% omegaOMA =  [10.0667;
%    31.4607;
%    49.1633;
%    63.2812;
%    72.8966];
% FDD
% omegaOMA =  [10.3544;
%    31.6000;
%    49.6243;
%    63.5835;
%    72.9024];
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
J=sum((omegas-omegaOMA).^2);
end