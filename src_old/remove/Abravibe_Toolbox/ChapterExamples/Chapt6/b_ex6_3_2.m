% b_ex6_3_2         Example 6.3.2 modal mass and stiffness
%
% Values are echoed to the screen by this example.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
clc
close all

% Following is repeated from a_ex6_3_1.m
M=eye(2);
K=[250 -150;-150 250];
A=inv(M)*K;
[V,D]=eig(A)

% Now to example 6.3.2:
% Modal mass:
ModalMass=V'*M*V

% Modal Stiffness:
ModalStiffness=V'*K*V

% Rescale mode shapes to unity modal stiffness:
Vnew(:,1)=V(:,1)/10;
Vnew(:,2)=V(:,2)/20

% Check:
NewModalStiffness=Vnew'*K*Vnew