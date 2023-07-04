% h_FRF_summary synthesizing FRFs from modal parameters
%
% Note that this example uses the system with prop. damping from Example
% 6.3.4, with the mode shapes scaled to unity modal mass. 

% This example shows all calculations to obtain FRFs, starting from
% the undamped system to obtain mode shapes, then changing to complex poles
% using the mode shape orthogonality. After that, the modal model is used
% to compute the frequency responses. The partial fraction expansion in
% this example is made using the MATLAB RESIDUE command.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
close all
clc

% Set up system
M=eye(2);
k1=100;
k2=150;
k3=100;
K=[k1+k2 -k2;-k2 k2+k3];
a=2/15;
b=1/1500;
C=a*M+b*K;

% Basic undamped calculations
A=inv(M)*K
[V,D]=eig(A)
wn=sqrt(diag(D))
Lambda=sqrt(-diag(D))       % These are the poles of the undamped system

% Weighted orthogonality criteria gives us the modal mass, stiffness and
% damping
Mr=diag(V'*M*V)
Kr=diag(V'*K*V)
Cr=diag(V'*C*V)

% Damped poles from modal coordinates
zeta=Cr./(2*sqrt(Mr.*Kr))
poles=-zeta.*wn+j*wn.*sqrt(1-zeta.^2)

%==================================================
% Residues
% Define the two points (p and q in the book)
p1=1;                   % Point 1
p2=1;                   % Point 2
% Compute the numerator and denominator of the transfer function, for the 
% first mode; this essentially comes from one line (the first) from the 
% Laplace transform of Eq. (6.62)
B=V(p1,1)*V(p2,1)/Mr(1)         % Numerator
A=[Mr(1) Cr(1) Kr(1)]           % Denominator
[R1,P1]=residue(B,A)            % Partial fraction expansion
% and the same thing for the second mode
B=V(p1,2)*V(p2,2)/Mr(2)
A=[Mr(2) Cr(2) Kr(2)]
[R2,P2]=residue(B,A)
% Form the FRF by summing over residues/poles
f=(0:.001:5)';
jw=j*2*pi*f;
% First driving point, H11
H(:,1)=R1(1)./(jw-P1(1))+R1(2)./(jw-P1(2))+R2(1)./(jw-P2(1))+R2(2)./(jw-P2(2));

% Now we redo the same procedure for FRF H12, i.e. force in DOF 2, response 
% in DOF 1
p1=1;                   % Point 1
p2=2;                   % Point 2
B=V(p1,1)*V(p2,1)/Mr(1)
A=[Mr(1) Cr(1) Kr(1)]
[R1,P1]=residue(B,A)
% and the second mode
B=V(p1,2)*V(p2,2)/Mr(2)
A=[Mr(2) Cr(2) Kr(2)]
[R2,P2]=residue(B,A)
H(:,2)=R1(1)./(jw-P1(1))+R1(2)./(jw-P1(2))+R2(1)./(jw-P2(1))+R2(2)./(jw-P2(2));

% Now we plot the two frequency responses on top of each other
subplot(2,1,1)
semilogy(f,abs(H))
legend('H_{11}','H_{12}')
ylabel('Dyn. Flexibility FRF [m/N]')
grid
subplot(2,1,2)
plot(f,angledeg(H))
ylabel('Phase, [Degrees]')
grid
xlabel('Frequency [Hz]')
