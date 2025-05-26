import numpy as np
import matplotlib



#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams



# Plot parameters
line = 2.5  # Line width
alp = 0.6  # alpha

rcParams.update({"figure.autolayout": True})
rcParams["axes.formatter.limits"] = [-3, 3]
rcParams["font.size"] = 12

# Loading dataset




vol_data= np.loadtxt("../data/volume_ints.dat")



start = 0



t = vol_data[start:, 0]
rho = vol_data[start:, 1]
source = vol_data[start:, 2]*0.8

# source_tot = np.cumsum(source)

# rho_diff = (- rho[0]  + rho) 


#BH_data = np.loadtxt("../data/stats_AH1.dat")


#BH_mass = BH_data[:,3]
# t2 = BH_data[:,0]





#net = (rho[0]-rho) - sum_flux

error = 100* (-rho[0] + rho)/(-rho[0])
fig, ax1 = plt.subplots(1, 1, figsize=(10,6))
width = 2.0

ax1.plot(t, error, linestyle='-', color='r', linewidth = width)
#ax1.plot(t, rho_diff, linestyle='-', color='r', linewidth = width)
#ax1.plot(t, source_tot, linestyle='-', color='b', linewidth = width)



# ax1.axvline(x=61.0, linestyle=':', color='k', linewidth = width)
# ax1.axhline(y=-rho[0], linestyle=':', color='k', linewidth = width)
#ax1.plot(t2, rho_diff2, linestyle='-', color='r')
#ax1.plot(t2, source_tot2, linestyle='--', color='r')




ax1.set_yscale("log")

ax1.set_xlabel("$t$", fontsize = 12)
# ax1.set_ylabel("$\\rho$", fontsize = 12)
#ax1.set_xlim(0,180)
#ax1.set_ylim(-0.3, 0.25)





#ax1.legend(["BH mass","$\\rho- \\rho_0$", "source"])



plt.savefig("../plots/error.png")