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

#%matplotlib qt

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
data_dir = pjoin(dirname(sio.__file__), r"C:\Users\Frede\OneDrive - Danmarks Tekniske Universitet\Github\Digital-Twin\src\python\data\newmark_jan_damp_simulation")
#data_dir = pjoin(dirname(sio.__file__), r"C:\Users\Christina\OneDrive - Danmarks Tekniske Universitet (1)\Github\Digital-Twin\src\data")
#data_dir = pjoin(dirname(sio.__file__), r"C:\Users\User\OneDrive\Dokumenter\GitHub\Digital-Twin\src\data")

#data_sim = pjoin(data_dir, 'data_sim_newmark_jan.mat')
Omega = []
Modes = []
Damp = []
for i in range(1000):
    data_sim = pjoin(data_dir, f'{i+1}_data_sim_newmark_jan_damp.mat')
    
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
    
    # ======== ANALYSIS ===========================================================
    # Run SSI
    br = 260
    SSIcov,Result = oma.SSIcovStaDiag(data, fs, br, ordmax=100)
    
    # Frequencies ccoming from the stability diagram
    # Anela
    #FreQ = [1.65612, 5.02666, 7.90311, 10.1157, 11.5903]
    # Jan new damp
    FreQ = [1.74726, 5.18409, 8.11932, 10.2961, 11.6574]
    # Jan new without damp
    #FreQ = [1.7473, 5.18712, 8.12024, 10.2972, 11.6468]
    # Extract the modal properties
    Res_SSIcov = oma.SSIModEX(FreQ, Result)
    
    # =============================================================================
    # Make some plots
    # =============================================================================
    MS_SSIcov = Res_SSIcov['Mode Shapes'].real
    
    Omega = np.append(Omega,Res_SSIcov["Frequencies"])
    Modes.append(MS_SSIcov)
    Damp = np.append(Damp,Res_SSIcov["Damping"])
    print(i)
    
np.save("data\Modal_parameters\Omegas",Omega)
np.save("data\Modal_parameters\Modes",Modes)
np.save("data\Modal_parameters\Damp",Damp)
#np.save(f"data\Modal_parameters\{i+1}_stab",Result['Reduced Poles'])
# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))