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

dataset = float(input('Which data set (1 (data_5_2_1)? and 2 (data_nodamp)? '))
# To open a .txt file create a variable containing the path to the file
if dataset == 1:
    _file = r"C:\Users\Frede\Speciale\Digital-Twin\src_operational\data\experimental_data\data_5_2_1.txt" # Path to the txt file

elif dataset == 2:
    _file = r"C:\Users\Frede\Speciale\Digital-Twin\src_operational\data\experimental_data\data.txt" # Path to the txt file




#_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_6_1.txt" # Path to the txt file
#_file = r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Christina\OneDrive - Danmarks Tekniske Universitet (1)\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data = pd.read_csv(_file, header=0, delim_whitespace=True, index_col=False) 
data = data.to_numpy()

data_base = data[:,6]

# remove all other columns
data = data[:,[1,2,3,4,5]]

# Swap the sensors because the sensor 1 is at the top
data = np.flip(data,1)

data_all_sens = data

data = data[:,[1,2,4]]

data = data #- data_base[:,None]

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
#First 2 modes
sos = signal.butter(75, 6.5/fs*2, btype='lowpass', output='sos')
filtdata_0_cut = signal.sosfiltfilt(sos, data.T) # filtered data

# Last 3 modes
sos = signal.butter(75, [6.5/fs*2, 14/fs*2], btype='bandpass', output='sos')
filtdata_cut_end = signal.sosfiltfilt(sos, data.T) # filtered data


# Filter
#First 1 modes
sos = signal.butter(75, (1.6570+5.01123)/2/fs*2, btype='lowpass', output='sos')
filtdata_mode1 = signal.sosfiltfilt(sos, data.T) # filtered data

# 2 mode
sos = signal.butter(75, [(1.6570+5.0170)/2/fs*2, (5.0170+7.8989)/2/fs*2], btype='bandpass', output='sos')
filtdata_mode2 = signal.sosfiltfilt(sos, data.T) # filtered data

# 3 mode
sos = signal.butter(75, [(5.01123+7.8989)/2/fs*2, (7.8989+10.1123)/2/fs*2], btype='bandpass', output='sos')
filtdata_mode3 = signal.sosfiltfilt(sos, data.T) # filtered data

# 4 mode
sos = signal.butter(75, [(7.8989+10.1123)/2/fs*2, (10.1123+11.5961)/2/fs*2], btype='bandpass', output='sos')
filtdata_mode4 = signal.sosfiltfilt(sos, data.T) # filtered data

# 5 mode
sos = signal.butter(75, [(10.1123+11.5961)/2/fs*2, 14/fs*2], btype='bandpass', output='sos')
filtdata_mode5 = signal.sosfiltfilt(sos, data.T) # filtered data


# Filter
sos = signal.butter(75, 14/fs*2, btype='lowpass', output='sos')
filtdata_20 = signal.sosfiltfilt(sos, data.T) # filtered data

sos = signal.butter(75, 6.5/fs*2, btype='lowpass', output='sos')
filtdata_9 = signal.sosfiltfilt(sos, data.T) # filtered data

# saving file with all sensors which are filtered, detrended and decimated
# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
fs=100
data_all_sens = signal.detrend(data_all_sens, axis=0) # Trend rmoval
q = 2 # Decimation factor
data_all_sens = signal.decimate(data_all_sens,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency
# Filter
sos = signal.butter(75, 14/fs*2, btype='lowpass', output='sos')
filtdata_all_sens = signal.sosfiltfilt(sos, data_all_sens.T) # filtered data

sos = signal.butter(75, 6.5/fs*2, btype='lowpass', output='sos')
filtdata_all_sens_2_first_modes = signal.sosfiltfilt(sos, data.T) # filtered data

if dataset == 1:
    np.save("Filtered_data_3_sensors\data_filt_3_sensors_high",data)
    np.save("Filtered_data_3_sensors\data_filt_0_cut_3_sensors_high",filtdata_0_cut)
    np.save("Filtered_data_3_sensors\data_filt_cut_end_3_sensors_high",filtdata_cut_end)
    np.save("Filtered_data_3_sensors\data_filt_mode1_3_sensors_high",filtdata_mode1)
    np.save("Filtered_data_3_sensors\data_filt_mode2_3_sensors_high",filtdata_mode2)
    np.save("Filtered_data_3_sensors\data_filt_mode3_3_sensors_high",filtdata_mode3)
    np.save("Filtered_data_3_sensors\data_filt_mode4_3_sensors_high",filtdata_mode4)
    np.save("Filtered_data_3_sensors\data_filt_mode5_3_sensors_high",filtdata_mode5)
    np.save("Filtered_data_3_sensors\data_filt_all_3_sensors_high",filtdata_20)
    np.save("Filtered_data_3_sensors\data_filt_first_modes_3_sensors_high",filtdata_9)
    np.save("Filtered_data_3_sensors\data_filt_all_sensors_high",filtdata_all_sens)    
    np.save("Filtered_data_3_sensors\data_filt_first_modes_all_sensorss_high",filtdata_all_sens_2_first_modes)  

elif dataset == 2:
    np.save("Filtered_data_3_sensors\data_filt_3_sensors_no_damp",data)
    np.save("Filtered_data_3_sensors\data_filt_0_cut_3_sensors_no_damp",filtdata_0_cut)
    np.save("Filtered_data_3_sensors\data_filt_cut_end_3_sensors_no_damp",filtdata_cut_end)
    np.save("Filtered_data_3_sensors\data_filt_mode1_3_sensors_no_damp",filtdata_mode1)
    np.save("Filtered_data_3_sensors\data_filt_mode2_3_sensors_no_damp",filtdata_mode2)
    np.save("Filtered_data_3_sensors\data_filt_mode3_3_sensors_no_damp",filtdata_mode3)
    np.save("Filtered_data_3_sensors\data_filt_mode4_3_sensors_no_damp",filtdata_mode4)
    np.save("Filtered_data_3_sensors\data_filt_mode5_3_sensors_no_damp",filtdata_mode5)
    np.save("Filtered_data_3_sensors\data_filt_all_3_sensors_no_damp",filtdata_20)
    np.save("Filtered_data_3_sensors\data_filt_first_modes_3_sensors_no_damp",filtdata_9)
    np.save("Filtered_data_3_sensors\data_filt_all_sensors_no_damp",filtdata_all_sens) 
    np.save("Filtered_data_3_sensors\data_filt_first_modes_all_sensors_no_damp",filtdata_all_sens_2_first_modes)  
    
#%% 
# plot plot plot
f, Pxx = signal.welch(data[:,0], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_low, Pxx_low = signal.welch(filtdata_0_cut[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_high, Pxx_high = signal.welch(filtdata_cut_end[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')

plt.figure(figsize=(10, 5))
plt.plot(f,20*np.log10(Pxx),label='Simulated response')
plt.plot(f_low,20*np.log10(Pxx_low),label='Lower frequency response')
plt.plot(f_high,20*np.log10(Pxx_high),label='Higher frequency response')
plt.ylabel('PSD [dB rel. to unit^2]')
plt.xlabel('Frequency [Hz]')
plt.legend()
plt.grid(True)

#%% 
f, Pxx = signal.welch(data[:,0], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_1, Pxx_1 = signal.welch(filtdata_mode1[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_2, Pxx_2 = signal.welch(filtdata_mode2[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_3, Pxx_3 = signal.welch(filtdata_mode3[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_4, Pxx_4 = signal.welch(filtdata_mode4[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')
f_5, Pxx_5 = signal.welch(filtdata_mode5[0,:], fs, nperseg=2**13, detrend='linear', return_onesided=True, scaling='spectrum')

plt.figure(figsize=(10, 5))
plt.plot(f,20*np.log10(Pxx),label='Structural response')
plt.plot(f_1,20*np.log10(Pxx_1),label='Frequency response for first mode')
plt.plot(f_2,20*np.log10(Pxx_2),label='Frequency response for second mode')
plt.plot(f_3,20*np.log10(Pxx_3),label='Frequency response for third mode')
plt.plot(f_4,20*np.log10(Pxx_4),label='Frequency response for fouth mode')
plt.plot(f_5,20*np.log10(Pxx_5),label='Frequency response for five mode')
plt.ylabel('PSD [dB rel. to unit^2]')
plt.xlabel('Frequency [Hz]')
plt.legend()
plt.grid(True)

Pxx = 20*np.log10(Pxx)
Pxx_1 = 20*np.log10(Pxx_1)
Pxx_2 = 20*np.log10(Pxx_2)
Pxx_3 = 20*np.log10(Pxx_3)
Pxx_4 = 20*np.log10(Pxx_4)
Pxx_5 = 20*np.log10(Pxx_5)

if dataset == 1:
    np.save("Filtered_data_3_sensors\data_filt_f_3_sensors_high",f)
    np.save("Filtered_data_3_sensors\data_filt_f_1_3_sensors_high",f_1)
    np.save("Filtered_data_3_sensors\data_filt_f_2_3_sensors_high",f_2)
    np.save("Filtered_data_3_sensors\data_filt_f_3_3_sensors_high",f_3)
    np.save("Filtered_data_3_sensors\data_filt_f_4_3_sensors_high",f_4)
    np.save("Filtered_data_3_sensors\data_filt_f_5_3_sensors_high",f_5)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_3_sensors_high",Pxx)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_1_3_sensors_high",Pxx_1)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_2_3_sensors_high",Pxx_2)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_3_3_sensors_high",Pxx_3)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_4_3_sensors_high",Pxx_4)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_5_3_sensors_high",Pxx_5)
elif dataset == 2:
    np.save("Filtered_data_3_sensors\data_filt_f_3_sensors_no_damp",f)
    np.save("Filtered_data_3_sensors\data_filt_f_1_3_sensors_no_damp",f_1)
    np.save("Filtered_data_3_sensors\data_filt_f_2_3_sensors_no_damp",f_2)
    np.save("Filtered_data_3_sensors\data_filt_f_3_3_sensors_no_damp",f_3)
    np.save("Filtered_data_3_sensors\data_filt_f_4_3_sensors_no_damp",f_4)
    np.save("Filtered_data_3_sensors\data_filt_f_5_3_sensors_no_damp",f_5)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_3_sensors_no_damp",Pxx)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_1_3_sensors_no_damp",Pxx_1)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_2_3_sensors_no_damp",Pxx_2)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_3_3_sensors_no_damp",Pxx_3)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_4_3_sensors_no_damp",Pxx_4)
    np.save("Filtered_data_3_sensors\data_filt_Pxx_5_3_sensors_no_damp",Pxx_5)
    