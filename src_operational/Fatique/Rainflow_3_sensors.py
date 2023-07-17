# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:41:41 2023

@author: Frede
"""

import numpy as np
import pandas as pd
import rainflow 

dataset = float(input('Which data set (1 (data_5_2_1)? and 2 (data_nodamp)? '))
# To open a .txt file create a variable containing the path to the file
if dataset == 1:
    data_me = np.load('C:/Users/Frede/OneDrive - Danmarks Tekniske Universitet/Speciale/Digital-Twin/src_operational/Virtual_sensing/Experimental/Filtered_data_3_sensors/data_filt_all_3_sensors_high.npy')
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src_operational\Virtual_sensing\Experimental\Filtered_data_3_sensors\data_predicted_3_sensors_high.txt" # Path to the txt file 
elif dataset == 2:
    data_me = np.load('C:/Users/Frede/OneDrive - Danmarks Tekniske Universitet/Speciale/Digital-Twin/src_operational/Virtual_sensing/Experimental/Filtered_data_3_sensors/data_filt_all_3_sensors_no_damp.npy')
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src_operational\Virtual_sensing\Experimental\Filtered_data_3_sensors\data_predicted_3_sensors_nodamp.txt" # Path to the txt file     
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data_pr = pd.read_csv(_file, header=0, delim_whitespace=True, index_col=False) 
data_pr = data_pr.to_numpy()
data_pr=data_pr[:,0]

# transpose the data
data_me=data_me.T
data_me=data_me[:,0]

rain_pr_binsize = rainflow.count_cycles(data_pr,binsize=1)
rain_me_binsize = rainflow.count_cycles(data_me,binsize=1)



if dataset == 1:
    np.savetxt('.\Rainflow_data_3_sensors\pr_rainflow_3_sensors_high.csv', rain_pr_binsize, delimiter=';')
    np.savetxt(".\Rainflow_data_3_sensors\me_rainflow_3_sensors_high.csv", rain_me_binsize, delimiter=';')
elif dataset == 2:
    np.savetxt('.\Rainflow_data_3_sensors\pr_rainflow_3_sensors_nodamp.csv', rain_pr_binsize, delimiter=';')
    np.savetxt('.\Rainflow_data_3_sensors\me_rainflow_3_sensors_nodamp.csv', rain_me_binsize, delimiter=';')