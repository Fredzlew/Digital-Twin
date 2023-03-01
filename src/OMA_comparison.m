clc;
clear all;
close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'))
% Loading the numerical model ad the simulated data
num = load('modelprop.mat');
SSI = load('SSImodalsim.mat');
SSI_new = load('SSImodalsim_newmark.mat');
ERA = load('ERAmodalsim.mat');
ERA_new = load('ERAmodalsim_newmark.mat');
FDD = load('FDDmodalsim.mat');
FDD_new = load('FDDmodalsim_newmark.mat');

% numerical modal parameters
numFreq = num.fn;
numphi = num.U;

% SSI modal parameters
SSIFreq = SSI.SSIFreq;
SSIphi = SSI.phi_SSI;

% SSI for newmark modal parameters
SSIFreq_new = SSI_new.SSIFreq;
SSIphi_new = SSI_new.phi_SSI;

% ERA modal parameters
ERAFreq = ERA.ERAFreq;
ERAphi = ERA.phi_ERA;

% ERA for newmark modal parameters
ERAFreq_new = ERA_new.ERAFreq;
ERAphi_new = ERA_new.phi_ERA;

% FDD modal parameters
FDDFreq = FDD.fn';
FDDphi = FDD.phi_FDD;

% FDD for newmark modal parameters
FDDFreq_new = FDD_new.fn';
FDDphi_new = FDD_new.phi_FDD;

% Comparing the data
% first all frequencies
modes = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5'};
Freq = table(numFreq,SSIFreq,ERAFreq,FDDFreq,SSIFreq_new,ERAFreq_new,FDDFreq_new,...
    'RowNames',modes);
disp('The frequencies out from simulated data  :')
disp(Freq)

modestot = {'Mode 1';'Mode 2';'Mode 3';'Mode 4';'Mode 5';'Total'};
q = 1;
for i = 1:5
    SSIdiff1(i) = q - min(SSIFreq(i),numFreq(i))/max(SSIFreq(i),numFreq(i));
    ERAdiff1(i) = q - min(ERAFreq(i),numFreq(i))/max(ERAFreq(i),numFreq(i));
    FDDdiff1(i) = q - min(ERAFreq(i),numFreq(i))/max(ERAFreq(i),numFreq(i));
    SSIdiff_new1(i) = q - min(SSIFreq_new(i),numFreq(i))/max(SSIFreq_new(i),numFreq(i));
    ERAdiff_new1(i) = q - min(ERAFreq_new(i),numFreq(i))/max(ERAFreq_new(i),numFreq(i));
    FDDdiff_new1(i) = q - min(FDDFreq_new(i),numFreq(i))/max(FDDFreq_new(i),numFreq(i));
end
SSIdiff = [SSIdiff1';sum(SSIdiff1)];
ERAdiff = [ERAdiff1';sum(ERAdiff1)];
FDDdiff = [FDDdiff1';sum(FDDdiff1)];
SSIdiff_new = [SSIdiff_new1';sum(SSIdiff_new1)];
ERAdiff_new = [ERAdiff_new1';sum(ERAdiff_new1)];
FDDdiff_new = [FDDdiff_new1';sum(FDDdiff_new1)];

% Difference betweeen the frequencies out from the numerical model
diff = table(SSIdiff,ERAdiff,FDDdiff,SSIdiff_new,ERAdiff_new,FDDdiff_new,...
    'RowNames',modestot);

disp('Difference between numerical frequencies and OMA frequencies')
disp(diff)