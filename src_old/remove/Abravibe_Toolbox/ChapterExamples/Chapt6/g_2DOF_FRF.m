% g_2DOF_FRF        FRF from 2DOF system
%
% Calculate frequency response of 2DOF system, from the mass, stiffness and
% matrices and modal damping. This command uses the ABRAVIBE toolbox command
% MKZ2FRF.

% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.

clear
clc
close all

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------

% Following is repeated from a_ex6_3_1.m
M=eye(2);
K=[250 -150;-150 250];
z=[0.01 0.01];                 % Damping of each mode

% Create a suitable frequency axis for computation of the frequency
% response
f=(0:5/2049:5-1/2049)';
% Now calculate the FRF (dyn. flexibility) between force in DOF 1, and response in DOF 2
Hd=mkz2frf(f,M,K,z,1,2,'d');

% Plot frequency response
semilogy(f,abs(Hd),'LineWidth',LineWidth)
xlabel('Frequency [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Frequency Response [m/N]','FontName',FontName,'FontSize',FontSize)
% axis([0 .5 -1e-3 1e-3])
grid
set(gca,'FontName',FontName,'FontSize',FontSize)

