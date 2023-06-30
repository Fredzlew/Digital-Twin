# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 09:28:35 2023

@author: Frede
"""


# Import modules
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import signal
import matplotlib.pyplot as plt
import PyOMA as oma
from matplotlib.backend_bases import MouseButton
from os.path import dirname, join as pjoin
import scipy.io as sio



# ======== PRE-PROCESSING =====================================================

# To open a .txt file create a variable containing the path to the file
# ======== PRE-PROCESSING =====================================================

# To open a .txt file create a variable containing the path to the file

data_dir = pjoin(dirname(sio.__file__), r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Speciale\Digital-Twin\src\python\data\newmark_jan_damp_simulation")





# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data_sim = pjoin(data_dir, '1_data_sim_newmark_jan_damp.mat')
#data_sim = pjoin(data_dir, f'{i+1}_data_sim_newmark_jan_damp.mat')
# Loading simulated data instead:
data = sio.loadmat(data_sim)
data = np.transpose(data["dis_new"])


# Sampling frequency
fs = 1000 # [Hz] Sampling Frequency

# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
q = 2 # Decimation factor
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

np.save("data\sim_filtered\sim_data",data)

# Filter
#First 3 modes
sos1 = signal.butter(4, 6.5/fs*2, btype='lowpass', output='sos')
filtdata_0_cut = signal.sosfiltfilt(sos1, data.T) # filtered data

# Last 2 modes
sos2 = signal.butter(4, 6.5/fs*2, btype='highpass', output='sos')
filtdata_cut_end = signal.sosfiltfilt(sos2, data.T) # filtered data

np.save("data\sim_filtered\data_filt_0_cut_sim",filtdata_0_cut)
np.save("data\sim_filtered\data_filt_cut_end_sim",filtdata_cut_end)


sos3 = signal.butter(4, 20/fs*2, btype='lowpass', output='sos')
filtdata_20 = signal.sosfiltfilt(sos3, data.T) # filtered data

sos4 = signal.butter(4, 6.5/fs*2, btype='lowpass', output='sos')
filtdata_9 = signal.sosfiltfilt(sos4, data.T) # filtered data


np.save("data\sim_filtered\data_filtdata_20",filtdata_20)
np.save("data\sim_filtered\data_filtdata_9",filtdata_9)

