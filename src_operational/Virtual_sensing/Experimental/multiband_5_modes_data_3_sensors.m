%% Script that performs the multiband %%
clear;clc;close all

% Adding path to data
addpath(genpath('.\Filtered_data_3_sensors'))

promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);
if q == 1
    % Load data
    f = readNPY('.\Filtered_data_3_sensors\data_filt_f_3_sensors_high.npy');
    f_1 = readNPY('.\Filtered_data_3_sensors\data_filt_f_1_3_sensors_high.npy');
    f_2 = readNPY('.\Filtered_data_3_sensors\data_filt_f_2_3_sensors_high.npy');
    f_3 = readNPY('.\Filtered_data_3_sensors\data_filt_f_3_3_sensors_high.npy');
    f_4 = readNPY('.\Filtered_data_3_sensors\data_filt_f_4_3_sensors_high.npy');
    f_5 = readNPY('.\Filtered_data_3_sensors\data_filt_f_5_3_sensors_high.npy');
    Pxx = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_3_sensors_high.npy');
    Pxx_1 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_1_3_sensors_high.npy');
    Pxx_2 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_2_3_sensors_high.npy');
    Pxx_3 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_3_3_sensors_high.npy');
    Pxx_4 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_4_3_sensors_high.npy');
    Pxx_5 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_5_3_sensors_high.npy');
    T = array2table([num2cell(f(1:2295)),num2cell(f_1(1:2295)),num2cell(f_2(1:2295)),num2cell(f_3(1:2295)),num2cell(f_4(1:2295)),num2cell(f_5(1:2295)),num2cell(Pxx(1:2295)),num2cell(Pxx_1(1:2295)),num2cell(Pxx_2(1:2295)),num2cell(Pxx_3(1:2295)),num2cell(Pxx_4(1:2295)),num2cell(Pxx_5(1:2295))]);
    T.Properties.VariableNames(1:12) = {'f','f_1','f_2','f_3','f_4','f_5','Pxx','Pxx_1','Pxx_2','Pxx_3','Pxx_4','Pxx_5'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_multi5_data_3_sensors_highdamp.csv','Delimiter',';')
elseif q == 2
    % Load data
    f = readNPY('.\Filtered_data_3_sensors\data_filt_f_3_sensors_no_damp.npy');
    f_1 = readNPY('.\Filtered_data_3_sensors\data_filt_f_1_3_sensors_no_damp.npy');
    f_2 = readNPY('.\Filtered_data_3_sensors\data_filt_f_2_3_sensors_no_damp.npy');
    f_3 = readNPY('.\Filtered_data_3_sensors\data_filt_f_3_3_sensors_no_damp.npy');
    f_4 = readNPY('.\Filtered_data_3_sensors\data_filt_f_4_3_sensors_no_damp.npy');
    f_5 = readNPY('.\Filtered_data_3_sensors\data_filt_f_5_3_sensors_no_damp.npy');
    Pxx = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_3_sensors_no_damp.npy');
    Pxx_1 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_1_3_sensors_no_damp.npy');
    Pxx_2 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_2_3_sensors_no_damp.npy');
    Pxx_3 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_3_3_sensors_no_damp.npy');
    Pxx_4 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_4_3_sensors_no_damp.npy');
    Pxx_5 = readNPY('.\Filtered_data_3_sensors\data_filt_Pxx_5_3_sensors_no_damp.npy');
    T = array2table([num2cell(f(1:2295)),num2cell(f_1(1:2295)),num2cell(f_2(1:2295)),num2cell(f_3(1:2295)),num2cell(f_4(1:2295)),num2cell(f_5(1:2295)),num2cell(Pxx(1:2295)),num2cell(Pxx_1(1:2295)),num2cell(Pxx_2(1:2295)),num2cell(Pxx_3(1:2295)),num2cell(Pxx_4(1:2295)),num2cell(Pxx_5(1:2295))]);
    T.Properties.VariableNames(1:12) = {'f','f_1','f_2','f_3','f_4','f_5','Pxx','Pxx_1','Pxx_2','Pxx_3','Pxx_4','Pxx_5'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_multi5_data_3_sensors_nodamp.csv','Delimiter',';')
end



