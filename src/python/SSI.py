
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

data_sim = pjoin(data_dir, 'data_sim_newmark_jan.mat')

# Loading simulated data instead:
data = sio.loadmat(data_sim)
data = np.transpose(data["dis_new"])


# Sampling frequency
fs = 1000 # [Hz] Sampling Frequency

# Using SciPy's signal module we can pre-process our data e.g. performing
# decimation, trend removal and filtering. 
# Detrend and decimate
data = signal.detrend(data, axis=0) # Trend rmoval
q = 10 # Decimation factor
data = signal.decimate(data,  q, ftype='fir', axis=0) # Decimation
fs = fs/q # [Hz] Decimated sampling frequency

# ======== ANALYSIS ===========================================================
# Run SSI
br = 130
SSIcov,Result = oma.SSIcovStaDiag(data, fs, br)

# Frequencies ccoming from the stability diagram
# Anela
#FreQ = [1.65612, 5.02666, 7.90311, 10.1157, 11.5903]
# Jan IRF
#FreQ = [1.74735, 5.18789, 8.1239, 10.3009, 11.6527]
# Jan new
FreQ = [1.7473, 5.18712, 8.12024, 10.2972, 11.6468]
# Extract the modal properties
Res_SSIcov = oma.SSIModEX(FreQ, Result)

# =============================================================================
# Make some plots
# =============================================================================
MS_SSIcov = Res_SSIcov['Mode Shapes'].real


np.save("omega",Res_SSIcov["Frequencies"])
np.save("phi",MS_SSIcov)
np.save("damp",Res_SSIcov["Damping"])