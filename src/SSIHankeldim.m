nm = 5; %Number of modes
nrows=50*nm/0.5;     %more than 20 * number of modes
ncols=4/5*10000;    %more than 2/3 of No. of data
options = optimset('Display','iter','PlotFcns',@optimplotfval);
[nrows,ncols] = fminsearch(@costfunSSIHankel,[rndeven(nrows),rndeven(ncols)],options);




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
nm = 5; %number of modes
output=f; % Displacements
ncols=4/5*length(f); % More than 2/3*number of samples
nrows=50*nm; % More than 20*number of sensors
cut=2*nm;  % cut=4 -> 2 modes, cut=10 -> 5 modes
[Result]=SSID(output,fs,ncols,nrows,cut);    %SSI
disp(Result.Parameters.NaFreq)