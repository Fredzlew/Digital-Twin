% a_ex6_3_1     Example 6.3.1 eigenvalues and eigenvectors
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

% Define mechanical system with two DOFs. See Example 6.3.1, p 124 in Brandt.
M=eye(2);
K=[250 -150;-150 250];

% Compute the matrix A
A=M\K;

% Compute the eigenvalues and eigenvectors
[V,D]=eig(A)
% Positive poles: -poles are also poles of this system, as all poles come
% in complex conjugate pairs
poles=sqrt(diag(-D))        % Corrected 2011-10-03, -D to obtain imaginary poles

% As you see, the eigenvectors are automatically scaled to unity length by
% MATLAB.
V