# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 10:18:26 2023

@author: Frede
"""

# import the builtin time module
import time
def cal_fact(n):
    if n == 1:
        return n
    else:
        return n * cal_fact(n-1)
# Grab Currrent Time Before Running the Code
start = time.time()
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

%matplotlib qt

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

# Filter
#First 3 modes
_b, _a = signal.butter(4, (0.1,6), fs=fs, btype='bandpass')
filtdata_0_cut = signal.filtfilt(_b, _a, data,axis=0) # filtered data

# Last 2 modes
_b, _a = signal.butter(4, (6,12), fs=fs, btype='bandpass')
filtdata_cut_end = signal.filtfilt(_b, _a, data,axis=0) # filtered data

filtdata = filtdata_0_cut + filtdata_cut_end

# ======== ANALYSIS ===========================================================
# Run SSI
br = 350
#SSIcov,Result = oma.SSIcovStaDiag(data, fs, br, ordmax=100)
SSIcov,Result = oma.SSIcovStaDiag(filtdata, fs, br, ordmax=100)

# Frequencies ccoming from the stability diagram

# Anela data_5_2_1
FreQ = [1.65592, 5.02282, 7.89698, 10.1326, 11.5861]



# Extract the modal properties
Res_SSIcov = oma.SSIModEX(FreQ, Result)

# =============================================================================
# Make some plots
# =============================================================================
MS_SSIcov = Res_SSIcov['Mode Shapes'].real

# Saving the parameters

np.save("data\Modal_parameters_anela\SSIfreq_nodamp_split",Res_SSIcov["Frequencies"])
np.save("data\Modal_parameters_anela\SSIphi_nodamp_split",MS_SSIcov)
np.save("data\Modal_parameters_anela\SSIdamp_nodamp_split",Res_SSIcov["Damping"])
np.save("data\Modal_parameters_anela\SSIstab_nodamp_split",Result['Reduced Poles'])
np.save("data\Modal_parameters_anela\data_nodamp_filtdata_split",filtdata)

# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))