% parameters
clc; clear; close all;
addpath(genpath('.\data'))
% Loading modal parameters from OMA 
promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
xx = input(promptt);
if xx == 1
    data = readNPY('.\data\experimental_data\Data_high.npy');
    data_de_dec = readNPY('.\data\experimental_data\Data_de_dec_high.npy');
    filtdata = readNPY('.\data\experimental_data\Filtdata_high.npy')';
    Time = readNPY('.\data\experimental_data\Time_high.npy');
    Time_dec = readNPY('.\data\experimental_data\Time_dec_high.npy');
elseif xx == 2
    data = readNPY('.\data\experimental_data\Data_nodamp.npy');
    data_de_dec = readNPY('.\data\experimental_data\Data_de_dec_nodamp.npy');
    filtdata = readNPY('.\data\experimental_data\Filtdata_nodamp.npy')';
    Time = readNPY('.\data\experimental_data\Time_nodamp.npy');
    Time_dec = readNPY('.\data\experimental_data\Time_dec_nodamp.npy');
end

% save data to latex
if xx == 1
    T1 = array2table([num2cell(Time(1:1000)),num2cell(data(1:1000,5))]);
    Table2 = array2table([num2cell(Time_dec(1:1000)),num2cell(data_de_dec(1:1000,5)),num2cell(filtdata(1:1000,5))]);
    T1.Properties.VariableNames(1:2) = {'Time','data5'};
    Table2.Properties.VariableNames(1:3) = {'Time_dec','data_de_dec5','filtdata5'};
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap3_data_high.csv','Delimiter',';')
    writetable(Table2,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap3_data_de_filt_high.csv','Delimiter',';')
elseif xx == 2
    T1 = array2table([num2cell(Time(1:1000)),num2cell(data(1:1000,5))]);
    Table2 = array2table([num2cell(Time_dec(1:1000)),num2cell(data_de_dec(1:1000,5)),num2cell(filtdata(1:1000,5))]);
    T1.Properties.VariableNames(1:2) = {'Time','data5'};
    Table2.Properties.VariableNames(1:3) = {'Time_dec','data_de_dec5','filtdata5'};
    writetable(T1,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap3_data_nodamp.csv','Delimiter',';')
    writetable(Table2,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap3_data_de_filt_nodamp.csv','Delimiter',';')
end
%% plot
plot(Time(1:1000),data(1:1000,5))