% i_3DOF_MKz        Create data from 3DOF system with M, K and z
%
% This example uses the simulation routines to compute FRFs and time data
% from a 3DOF model with known M and K matrices, and modal damping z

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
FontSize=11;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};


% Build M, K, and z
M=[1 0 0;0 2 0;0 0 1];
K=1e6*[2 -1 0;
    -1 2 -1;
    0 -1 2];
z=.01*[1 1 1];

% Build a modal model, with poles and mode shapes from M, K, and z:
[p,V]=mkz2modal(M,K,z);
[fr,zr]=poles2fz(p);

% Calculate and plot FRF (acceleration) between force and accel in DOF 1
fs=1000;
N=2*1024;
f=(0:fs/(8*N):fs/2)';
H=mkz2frf(f,M,K,z,1,1,'a');
semilogy(f,abs(H))
xlabel('Frequency, [Hz]','FontName',FontName,'FontSize',FontSize)
ylabel('Accelerance, [(m/s^2)/N]','FontName',FontName,'FontSize',FontSize)

% Generate input noise and compute output noise
N0=50;
x=sqrt(N0)*randn(500*N,1);
tic
y=timefresp(x,fs,p,V,1,1,'a');
T=toc;
Gyyt=N0/(fs/2)*abs(H).^2;

% Save the data in the data directory
Text='These data were created by i_3DOF_MKz.m in Chapter 6 directory';
save -mat ..\Data\3dofmkz.mat f Gyyt x y fs N Text M K z
fprintf('Data are saved in file ..\\Data\\3dofmkz\n')

%=========================================================================
% Here is a brief discussion of the accuracy of the timefresp command, and
% how to treat it:
% For most purposes, the approach above is accurate enough, the accuracy
% with 10 times oversampling timefresp uses is within 4% in the frequency
% range 0 to 0.4 fs.
% If you want 'absolute' accuracy, what you can do is to compute the true
% FRF of the system that the digital filters in timefresp define, rather
% than the FRF of the M, C, K system, which timefresp only manages to
% approximate. The actual FRF used by timefresp can be computed using the
% 'hidden' syntax
% >> [y,Htrue,ftrue]=timefresp(...);
% This produces the actual FRF of the digital filters used by the
% algorithm. See example b_SISORandnEx.m in the folder for chapter 13.