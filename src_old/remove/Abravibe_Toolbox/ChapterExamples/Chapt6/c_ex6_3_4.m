% c_ex6_3_4         Example 6.3.4 Proportional damping
%
% This example uses the mck2modal command from the ABRAVIBE toolbox to
% compute the modal parameters in a case of a system with proportional
% damping. See inside the file mck2modal for details!

% Ref: See Example 6.3.4 on page 128 ff in Brandt.

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
C=2/15*M+1/1500*K;

% Compute poles and mode shapes:
[p,V]=mck2modal(M,C,K)

% To read the poles in natural frequency and damping, use:
[fr,zr]=poles2fz(p)