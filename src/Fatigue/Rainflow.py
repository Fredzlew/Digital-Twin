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
    _file = r"..\Virtual_sensing\Experimental\Filtered_data\data_predicted_high.txt" # Path to the txt file 
elif dataset == 2:
    _file = r"..\Virtual_sensing\Experimental\Filtered_data\data_predicted_nodamp.txt" # Path to the txt file     
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data = pd.read_csv(_file, header=0, delim_whitespace=True, index_col=False) 
data = data.to_numpy()
data_pr=data[:,0]
data_me=data[:,1]


rain_pr_binsize = rainflow.count_cycles(data_pr,binsize=1)
rain_me_binsize = rainflow.count_cycles(data_me,binsize=1)



if dataset == 1:
    np.savetxt('.\Rainflow_data\pr_rainflow_high.csv', rain_pr_binsize, delimiter=';')
    np.savetxt(".\Rainflow_data\me_rainflow_high.csv", rain_me_binsize, delimiter=';')
elif dataset == 2:
    np.savetxt('.\Rainflow_data\pr_rainflow_nodamp.csv', rain_pr_binsize, delimiter=';')
    np.savetxt('.\Rainflow_data\me_rainflow_nodamp.csv', rain_me_binsize, delimiter=';')