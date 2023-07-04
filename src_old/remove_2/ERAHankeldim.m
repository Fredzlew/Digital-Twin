nm = 5; %Number of modes
nrows=50*(2*nm/5)+1;     %more than 20 * number of modes
ncols=4/5*10000-nrows-3;    %more than 2/3 of No. of data
options = optimset('Display','iter','PlotFcns',@optimplotfval);
[nrows,ncols] = fminsearch(@costfunERAHankel,[nrows,ncols],options);




% Choose of data
prompt = "Use ERA for measured or simulated data (1=measured, 2=simulated)? ";
ERAdata = input(prompt);
if ERAdata == 1
    % Measurements
    data = readmatrix('data_1_2_1.txt')'; % Loading displacement data
    fss = data(2:6,1:10000)/1000; % Converting mm to m
    f = [fss(5,:);fss(4,:);fss(3,:);fss(2,:);fss(1,:)]; % Swap columns due to sensor
elseif ERAdata == 2
    % Simulated data
    data_sim = load('data_sim.mat');
    f = data_sim.dis(:,1:10000);
end
fs=100; % Sampling frequency (1/dt)
Y=f; %Displacements
inputs=1;     
cut=2*nm;        %Identify 5 modes
shift=10;      %Adjust EMAC sensitivity
EMAC_option=1; %EMAC is calculated only from observability matrix

[Result]=ERA(Y,fs,ncols,nrows,inputs,cut,shift,EMAC_option);  %ERA
disp(Result.Parameters.NaFreq)