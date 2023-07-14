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
    data_me = np.load('C:/Users/Frede/OneDrive - Danmarks Tekniske Universitet/Speciale/Digital-Twin/src/Virtual_sensing/Experimental/Filtered_data/data_filt_all_high.npy')
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src\Virtual_sensing\Experimental\Filtered_data\data_predicted_high.txt" # Path to the txt file 
elif dataset == 2:
    data_me = np.load('C:/Users/Frede/OneDrive - Danmarks Tekniske Universitet/Speciale/Digital-Twin/src/Virtual_sensing/Experimental/Filtered_data/data_filt_all_no_damp.npy')
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src\Virtual_sensing\Experimental\Filtered_data\data_predicted_nodamp.txt" # Path to the txt file     
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data_pr = pd.read_csv(_file, header=0, delim_whitespace=True, index_col=False) 
data_pr = data_pr.to_numpy()
data_pr=data_pr[:,0]

# transpose the data
data_me=data_me.T
data_me=data_me[:,0]

#data=np.array(data).tolist()


rain_pr_nbin = rainflow.count_cycles(data_pr,nbins=30)
rain_me_nbin = rainflow.count_cycles(data_me,nbins=30)

rain_pr_binsize = rainflow.count_cycles(data_pr,binsize=1)
rain_me_binsize = rainflow.count_cycles(data_me,binsize=1)

rain_pr_ndig = rainflow.count_cycles(data_pr,ndigits=0)
rain_me_ndig = rainflow.count_cycles(data_me,ndigits=0)

#for rng, mean, count, i_start, i_end in rainflow.extract_cycles(data): 
    #print(1)