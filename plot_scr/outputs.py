

# Load the modules
import matplotlib
import yt
import time
import numpy as np
import math
import csv
import sys



matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams


start_time = time.time()
yt.enable_parallelism()

# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]
rcParams["font.size"] = 12

# CHANGE ME
# Loading dataset (Load from one folder up)
data_location = "../hdf5/DWp_*"  # Data file location


# Loading dataset
ts = yt.load(data_location)






my_storage = {}

for sto, ds in ts.piter(storage=my_storage):
  
   

   ad = ds.all_data()

   rho_max_val = ad.max("rho")
   
   
   
   array = [rho_max_val]
   
   sto.result = array
   sto.result_id = str(ds)
 

 

if yt.is_root():
     
     rho_max_val_dat = []

     for L in sorted(my_storage.items()):
       rho_max_val_dat.append(L[1][0])
     
       
       
     np.savetxt('../data/rho_max_val.out', rho_max_val_dat)
    


    


