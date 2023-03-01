clc;
clear all;
close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
% Loading the numerical model ad the simulated data
num = load('modelprop.mat');
SSI = load('SSImodalsim.mat');
ERA = load('ERAmodalsim.mat');
ERA_new = load('ERAmodalsim_newmark.mat');
FDD = load('FDDmodalsim.mat');

% numerical modal parameters
numFreq = num.fn;
numphi = num.U;

% SSI modal parameters
SSIFreq = SSI.SSIFreq;
SSIphi = SSI.phi_SSI;

% ERA modal parameters
ERAFreq = ERA.ERAFreq;
ERAphi = ERA.phi_ERA;

% ERA for newmark modal parameters
ERAFreq_new = ERA_new.ERAFreq;
ERAphi_new = ERA_new.phi_ERA;

% FDD modal parameters
FDDFreq = FDD.fn';
FDDphi = FDD.phi_FDD;

% Comparing the data
% first all frequencies
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
Freq = table(numFreq,SSIFreq,ERAFreq,FDDFreq,ERAFreqNew,...
    'RowNames',modes);
disp('The frequencies out from simulated data  :')
disp(Freq)

modestot = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5';'Total'};
SSIdiff1 = [100;100;100;100;100] - [min(SSIFreq(1),numFreq(1))/max(SSIFreq(1),numFreq(1));min(SSIFreq(2),numFreq(2))/max(SSIFreq(2),numFreq(2));min(SSIFreq(3),numFreq(3))/max(SSIFreq(3),numFreq(3));min(SSIFreq(4),numFreq(4))/max(SSIFreq(4),numFreq(4));min(SSIFreq(5),numFreq(5))/max(SSIFreq(5),numFreq(5))]*100;
ERAdiff1 = [100;100;100;100;100] - [min(ERAFreq(1),numFreq(1))/max(ERAFreq(1),numFreq(1));min(ERAFreq(2),numFreq(2))/max(ERAFreq(2),numFreq(2));min(ERAFreq(3),numFreq(3))/max(ERAFreq(3),numFreq(3));min(ERAFreq(4),numFreq(4))/max(ERAFreq(4),numFreq(4));min(ERAFreq(5),numFreq(5))/max(ERAFreq(5),numFreq(5))]*100;
FDDdiff1 = [100;100;100;100;100] - [min(FDDFreq(1),numFreq(1))/max(FDDFreq(1),numFreq(1));min(FDDFreq(2),numFreq(2))/max(FDDFreq(2),numFreq(2));min(FDDFreq(3),numFreq(3))/max(FDDFreq(3),numFreq(3));min(FDDFreq(4),numFreq(4))/max(FDDFreq(4),numFreq(4));min(FDDFreq(5),numFreq(5))/max(FDDFreq(5),numFreq(5))]*100;
ERAdiff_new1 = [100;100;100;100;100] - [min(ERAFreq_new(1),numFreq(1))/max(ERAFreq_new(1),numFreq(1));min(ERAFreq_new(2),numFreq(2))/max(ERAFreq_new(2),numFreq(2));min(ERAFreq_new(3),numFreq(3))/max(ERAFreq_new(3),numFreq(3));min(ERAFreq_new(4),numFreq(4))/max(ERAFreq_new(4),numFreq(4));min(ERAFreq_new(5),numFreq(5))/max(ERAFreq_new(5),numFreq(5))]*100;
SSIdiff = [SSIdiff1;sum(SSIdiff1)];
ERAdiff = [ERAdiff1;sum(ERAdiff1)];
ERAdiff_new = [ERAdiff_new1;sum(ERAdiff_new1)];
FDDdiff = [FDDdiff1;sum(FDDdiff1)];

% Difference betweeen the frequencies out from the numerical model
diff = table(SSIdiff,ERAdiff,FDDdiff,ERAdiff_new,...
    'RowNames',modestot);

disp('Difference between numerical frequencies and OMA frequencies')
disp(diff)