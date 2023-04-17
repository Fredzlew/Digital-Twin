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
# Create a figure and a set of subplots

%matplotlib qt

# ======== PRE-PROCESSING =====================================================
dataset = float(input('Which data set (1 (data_5_2_1) and 2 (data_5_6_1)? '))
# To open a .txt file create a variable containing the path to the file
if dataset == 1:
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_2_1.txt" # Path to the txt file

elif dataset == 2:
    _file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_5_6_1.txt" # Path to the txt file

#_file = r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Christina\OneDrive - Danmarks Tekniske Universitet (1)\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

data = pd.read_csv(_file, header=1, delim_whitespace=True, index_col=False) 
data = data.to_numpy()

# Removing the time column 
data = np.delete(data, 0, 1)
data = np.delete(data, 5, 1)

# Swap the sensors because the sensor 1 is at the top
data[:,[4,3,2,1,0]] = data[:,[0,1,2,3,4]]



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
_b, _a = signal.butter(4, (0.1,12), fs=fs, btype='bandpass')
filtdata = signal.filtfilt(_b, _a, data,axis=0) # filtered data

# ======== ANALYSIS ===========================================================
# Run FDD
#FDD, Result = oma.FDDsvp(data,  fs)
FDD, Result = oma.FDDsvp(filtdata,  fs)

# Finding frequency to plot the FDD
dt = 1/fs
N = data.shape[0]
f = np.transpose(np.linspace(0,N/2,num=115001))*dt

# Define list/array with the peaks identified from the plot
if dataset == 1:
    # Anela data_5_2_1
    FreQ = [1.66, 5.03, 7.90, 10.12, 11.59] # identified peaks

elif dataset == 2:
    # Anela data_5_6_1
    FreQ = [1.66, 5.030, 7.90, 10.12, 11.59] # identified peaks


# Extract the modal properties 
Res_EFDD = oma.EFDDmodEX(FreQ, Result, method='EFDD', MAClim=0.8)


MS_EFDD = Res_EFDD['Mode Shapes'].real

# Saving parameters
if dataset == 1:
    # Anela data_5_2_1
    np.save("data\Modal_parameters_anela\FDDomega_5_2_1",Res_EFDD["Frequencies"])
    np.save("data\Modal_parameters_anela\FDDphi_5_2_1",MS_EFDD )
    np.save("data\Modal_parameters_anela\FDDdamp_5_2_1",Res_EFDD["Damping"] )
    np.save("data\Modal_parameters_anela\FDDPSD_5_2_1",Result["Singular Values"])
    np.save("data\Modal_parameters_anela\FDDPSDfreq_5_2_1",f)

elif dataset == 2:
    # Anela data_5_6_1
    np.save("data\Modal_parameters_anela\FDDomega_5_6_1",Res_EFDD["Frequencies"])
    np.save("data\Modal_parameters_anela\FDDphi_5_6_1",MS_EFDD )
    np.save("data\Modal_parameters_anela\FDDdamp_5_6_1",Res_EFDD["Damping"] )
    np.save("data\Modal_parameters_anela\FDDPSD_5_6_1",Result["Singular Values"])
    np.save("data\Modal_parameters_anela\FDDPSDfreq_5_6_1",f)
    

# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))