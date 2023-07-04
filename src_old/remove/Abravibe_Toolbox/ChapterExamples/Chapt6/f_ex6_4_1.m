% f_ex6.4.1     Synthesizing FRFs from modal parameters
%
% Note that this example uses the system with prop. damping from Example
% 6.3.4, with the mode shapes scaled to unity modal mass. 

% This example shows all calculations to obtain FRFs, starting from
% the undamped system to obtain mode shapes, then changing to complex poles
% using the mode shape orthogonality. After that, the modal model is used
% to compute the frequency responses.

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

% Compute poles and mode shapes:
[p,V]=mck2modal(M,C,K)
% Rescale to unity modal mass, as mck2modal uses unity modal A
V=uma2umm(V,p);

% Modal scaling constants
Q1=1/(j*2*imag(p(1))*M(1,1))
Q2=1/(j*2*imag(p(2))*M(2,2))

% Form the FRF by summing over residues/poles. First we create a frequency
% axis
f=(0:.001:5)';          % In Hz
jw=j*2*pi*f;            % This is useful since denominators are in jw-lambda

% Now calculate the driving point FRF, H11.
% Residues
A111=Q1*V(1,1)*V(1,1)
A112=Q2*V(1,2)*V(1,2)
% and the FRF
H(:,1)=A111./(jw-p(1))+conj(A111)./(jw-conj(p(1)))+...
    A112./(jw-p(2))+conj(A112)./(jw-conj(p(2)));

% Now we redo the same procedure for FRF H12, i.e. force in DOF 2, response 
% in DOF 1
A121=Q1*V(2,1)*V(1,1)
A122=Q2*V(2,2)*V(1,2)
% and the FRF
H(:,2)=A121./(jw-p(1))+conj(A121)./(jw-conj(p(1)))+...
    A122./(jw-p(2))+conj(A122)./(jw-conj(p(2)));

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
