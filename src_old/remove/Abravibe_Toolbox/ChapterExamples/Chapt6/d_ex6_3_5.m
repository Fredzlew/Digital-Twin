% d_ex6_3_5         Example 6.3.5 Non-proportional damping
%
% This example uses the mck2modal command from the ABRAVIBE toolbox to
% compute the modal parameters in a case of a system with non-proportional
% damping. See inside the file mck2modal for details!

% Ref: See Example 6.3.5 on page 132 in Brandt.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
clc
close all

% Following is repeated from a_ex6_3_4.m
M=eye(2);
K=[250 -150;-150 250];
C=2/15*M+1/1500*K;
C(1,1)=C(1,1)+0.5;

%================================================================
% First, use the code from page 132 in the book:
A=[C M; M 0*M];
B=[K 0*M; 0*M -M];
[V,D]=eig(-A\B)
% Sort in descending order
[Dum,I]=sort(diag(abs(imag(D))));
p=diag(D(I,I))
V=V(:,I);
% Scale to unity Modal A
Ma=V.'*A*V;
for col = 1:length(V(1,:))
    Va(:,col)=V(:,col)/sqrt(Ma(col,col));
end
Ma_new=Va.'*A*Va;
fprintf('Mode shapes after scaling to unity modal A:\n')
Va

%================================================================
% Now compare with the result of the ABRAVIBE command MCK2MODAL:
% You may notice a sign shift in the second mode shape, also compared to the
% listings in the book. However, the modal scaling does not change with a 
% sign shift in a mode shape, since the mode shape always comes "multiplied
% twice"
fprintf('Using mck2modal:\n')
% Compute poles and mode shapes:
[p_cmd,V_cmd]=mck2modal(M,C,K)

% To read the poles in natural frequency and damping, use:
[fr,zr]=poles2fz(p_cmd)