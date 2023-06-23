function J=costfunERAHankel(dim)
nrows = dim(1);
ncols = dim(2);
% Choose of data
% prompt = 2;
% ERAdata = input(prompt);
% if ERAdata == 1
%     % Measurements
%     data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
%     fss = data(2:6,1:10000)/1000; % Converting mm to m
%     f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
% elseif ERAdata == 2
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
% end

filename = load('modelprop.mat'); % Loads mass and stiffness matrices
%Identify modal parameters using displacement with added uncertainty
%--------------------------------------------------------------------------
fs=100; % Sampling frequency (1/dt)
nm = 5; %Number of modes
Y=f; %Displacements
inputs=1;     
cut=2*nm;        %Identify 5 modes
shift=10;      %Adjust EMAC sensitivity
EMAC_option=1; %EMAC is calculated only from observability matrix

[Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA
% Find frequencies from ERA
omegaOMA = Result.Parameters.NaFreq;
% Find frequencies from modelprop
omegas = filename.fn;
% Cost function
J=sum((omegas-omegaOMA).^2);
disp(nrows)
disp(ncols)
end