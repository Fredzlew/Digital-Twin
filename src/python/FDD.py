
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

# To open a .txt file create a variable containing the path to the file
#_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
#_file = r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
#_file = r"C:\Users\Christina\OneDrive - Danmarks Tekniske Universitet (1)\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file
# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present

#data = pd.read_csv(_file, header=1, delim_whitespace=True, index_col=False) 
#data = data.to_numpy()

# Removing the time column 
#data = np.delete(data, 0, 1)
#data = np.delete(data, 5, 1)

# Swap the sensors because the sensor 1 is at the top
#data[:,[4,3,2,1,0]] = data[:,[0,1,2,3,4]]

#Finding the file
data_dir = pjoin(dirname(sio.__file__), r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data")
#data_dir = pjoin(dirname(sio.__file__), r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data")
#data_sim = pjoin(data_dir, 'data_sim_newmark_jan.mat')
data_sim = pjoin(data_dir, 'data_sim_newmark_jan_damp.mat')

# Loading simulated data instead:
data = sio.loadmat(data_sim)
data = np.transpose(data["dis_new"])



# Sampling frequency
fs = 1000 # [Hz] Sampling Frequency

# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
q = 5 # Decimation factor
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

# Filter
_b, _a = signal.butter(10, (0.3,6.5), fs=fs, btype='bandpass')
filtdata = signal.filtfilt(_b, _a, data,axis=0) # filtered data

# ======== ANALYSIS ===========================================================
# Run FDD
FDD, Result = oma.FDDsvp(data,  fs)
#FDD, Result = oma.FDDsvp(filtdata,  fs)

# Define list/array with the peaks identified from the plot
# Anela
#FreQ = [1.66012, 5.03053, 7.90999, 10.1107, 11.5936] # identified peaks
# Simuleret data
FreQ = [1.7508, 5.19052, 8.12081, 10.301, 11.6512] # identified peaks

# Extract the modal properties 
#Res_FDD = oma.FDDmodEX(FreQ, Result)
Res_EFDD = oma.EFDDmodEX(FreQ, Result, method='EFDD')
#Res_FSDD = oma.EFDDmodEX(FreQ, Result, method='FSDD', npmax = 35, MAClim=0.95, plot=True)

#MS_FDD = Res_FDD['Mode Shapes'].real
MS_EFDD = Res_EFDD['Mode Shapes'].real
#MS_FSDD = Res_FSDD['Mode Shapes'].real

np.save("FDDomega",Res_EFDD["Frequencies"])
np.save("FDDphi",MS_EFDD )
np.save("FDDdamp",Res_EFDD["Damping"] )