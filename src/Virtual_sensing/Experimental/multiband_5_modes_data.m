%% Script that performs the multiband %%
clear;clc;close all

% Adding path to data
addpath(genpath('.\Filtered_data'))

promptt = "High damping or no damping? (1 = High and 2 = no damp): ";
q = input(promptt);
if q == 1
    % Load data
    f = readNPY('.\Filtered_data\data_filt_f_high.npy');
    f_1 = readNPY('.\Filtered_data\data_filt_f_1_high.npy');
    f_2 = readNPY('.\Filtered_data\data_filt_f_2_high.npy');
    f_3 = readNPY('.\Filtered_data\data_filt_f_3_high.npy');
    f_4 = readNPY('.\Filtered_data\data_filt_f_4_high.npy');
    f_5 = readNPY('.\Filtered_data\data_filt_f_5_high.npy');
    Pxx = readNPY('.\Filtered_data\data_filt_Pxx_high.npy');
    Pxx_1 = readNPY('.\Filtered_data\data_filt_Pxx_1_high.npy');
    Pxx_2 = readNPY('.\Filtered_data\data_filt_Pxx_2_high.npy');
    Pxx_3 = readNPY('.\Filtered_data\data_filt_Pxx_3_high.npy');
    Pxx_4 = readNPY('.\Filtered_data\data_filt_Pxx_4_high.npy');
    Pxx_5 = readNPY('.\Filtered_data\data_filt_Pxx_5_high.npy');
    T = array2table([num2cell(f),num2cell(f_1),num2cell(f_2),num2cell(f_3),num2cell(f_4),num2cell(f_5),num2cell(Pxx),num2cell(Pxx_1),num2cell(Pxx_2),num2cell(Pxx_3),num2cell(Pxx_4),num2cell(Pxx_5)]);
    T.Properties.VariableNames(1:12) = {'f','f_1','f_2','f_3','f_4','f_5','Pxx','Pxx_1','Pxx_2','Pxx_3','Pxx_4','Pxx_5'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_multi5_data_highdamp.csv','Delimiter',';')
elseif q == 2
    % Load data
    f = readNPY('.\Filtered_data\data_filt_f_no_damp.npy');
    f_1 = readNPY('.\Filtered_data\data_filt_f_1_no_damp.npy');
    f_2 = readNPY('.\Filtered_data\data_filt_f_2_no_damp.npy');
    f_3 = readNPY('.\Filtered_data\data_filt_f_3_no_damp.npy');
    f_4 = readNPY('.\Filtered_data\data_filt_f_4_no_damp.npy');
    f_5 = readNPY('.\Filtered_data\data_filt_f_5_no_damp.npy');
    Pxx = readNPY('.\Filtered_data\data_filt_Pxx_no_damp.npy');
    Pxx_1 = readNPY('.\Filtered_data\data_filt_Pxx_1_no_damp.npy');
    Pxx_2 = readNPY('.\Filtered_data\data_filt_Pxx_2_no_damp.npy');
    Pxx_3 = readNPY('.\Filtered_data\data_filt_Pxx_3_no_damp.npy');
    Pxx_4 = readNPY('.\Filtered_data\data_filt_Pxx_4_no_damp.npy');
    Pxx_5 = readNPY('.\Filtered_data\data_filt_Pxx_5_no_damp.npy');
    T = array2table([num2cell(f),num2cell(f_1),num2cell(f_2),num2cell(f_3),num2cell(f_4),num2cell(f_5),num2cell(Pxx),num2cell(Pxx_1),num2cell(Pxx_2),num2cell(Pxx_3),num2cell(Pxx_4),num2cell(Pxx_5)]);
    T.Properties.VariableNames(1:12) = {'f','f_1','f_2','f_3','f_4','f_5','Pxx','Pxx_1','Pxx_2','Pxx_3','Pxx_4','Pxx_5'};
    writetable(T,'C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Kandidat\Data\Kap10_virtuel_sensing_multi5_data_nodamp.csv','Delimiter',';')
end



