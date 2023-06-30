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

_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data.txt" # Path to the txt file




#_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_6_1.txt" # Path to the txt file
#_file = r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Christina\OneDrive - Danmarks Tekniske Universitet (1)\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data = pd.read_csv(_file, header=0, delim_whitespace=True, index_col=False) 
data = data.to_numpy()

# remove all other columns
data = data[:,[1,2,3,4,5]]

# Swap the sensors because the sensor 1 is at the top
data = np.flip(data,1)


# Sampling frequency
fs = 100 # [Hz] Sampling Frequency

# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
q = 2 # Decimation factor
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

np.save("data\Modal_parameters_anela\data_filt_nodamp",data)

# Filter
#First 3 modes
sos = signal.butter(75, 9/fs*2, btype='lowpass', output='sos')
filtdata_0_cut = signal.sosfiltfilt(sos, data.T) # filtered data

# Last 2 modes
sos = signal.butter(75, 9/fs*2, btype='highpass', output='sos')
filtdata_cut_end = signal.sosfiltfilt(sos, data.T) # filtered data

np.save("data\Modal_parameters_anela\data_filt_0_cut",filtdata_0_cut)
np.save("data\Modal_parameters_anela\data_filt_cut_end",filtdata_cut_end)


sos = signal.butter(75, 20/fs*2, btype='lowpass', output='sos')
filtdata_20 = signal.sosfiltfilt(sos, data.T) # filtered data

sos = signal.butter(75, 9/fs*2, btype='lowpass', output='sos')
filtdata_9 = signal.sosfiltfilt(sos, data.T) # filtered data


np.save("data\Modal_parameters_anela\datafiltnodamp20",filtdata_20)
np.save("data\Modal_parameters_anela\data_filt_nodamp_6.5",filtdata_9)




