
# Import modules
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import signal
import matplotlib.pyplot as plt
import PyOMA as oma


# ======== PRE-PROCESSING =====================================================
# To open a .txt file create a variable containing the path to the file
_file = r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\data\Anela\data_1_2_1.txt" # Path to the txt file

# open the file with pandas and create a dataframe
# N.B. whatchout for header, separator and remove time column if present
data = pd.read_csv(_file, header=1, delim_whitespace=True, index_col=False) 
data = data.to_numpy()

# Removing the time column 
data = np.delete(data, 0, 1)
data = np.delete(data, 5, 1)

# Swap the sensors because the sensor 1 is at the top
data[:,[4,3,2,1,0]] = data[:,[0,1,2,3,4]]


# to retrieve the example data 
#data, (fex, FI_ex, xi_ex) = oma.Exdata()

# Sampling frequency
fs = 100 # [Hz] Sampling Frequency

"""
# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
q = 5 # Decimation factor
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

# Filter
_b, _a = signal.butter(12, (0.3,6.5), fs=fs, btype='bandpass')
filtdata = signal.filtfilt(_b, _a, data,axis=0) # filtered data

"""
# ======== ANALYSIS ===========================================================
# Run FDD
FDD = oma.FDDsvp(data,  fs)
# FDD = FDDsvp(filtdata,  fs)

# Define list/array with the peaks identified from the plot
FreQ = [1.64794921875,	5.029296875,	7.89794921875,	10.11962890625,	11.602783203125] # identified peaks

# Run SSI
br = 15
SSIcov= oma.SSIcovStaDiag(data, fs, br)

# Extract the modal properties
Res_SSIcov = oma.SSIModEX(FreQ, SSIcov[1])

# =============================================================================
# Make some plots
# =============================================================================
MS_SSIcov = Res_SSIcov['Mode Shapes'].real