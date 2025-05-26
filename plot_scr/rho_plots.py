import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np

import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
mpl.use("Agg")
import os
from unyt import unyt_array

F2 = 14 # Legend font size
line = 0.8 # Line width
alp = 1 # alpha
ticklabelsize = F2 


rc = mpl.rcParams # Font structure is called rc now
rc["text.usetex"] = True # Tex fonts
# rc["mathtext.fontset"] = "stix"
# rc["font.family"] = "serif"
# rc["font.sans-serif"] = "Helvetica"
# rc["font.serif"].insert(0,"cm") # Default font is computer modern for latex
rc["font.size"] = F2
rc["lines.linewidth"] = line
rc["axes.titlepad"] = 7.
rc["axes.axisbelow"] = False
rc["axes.linewidth"] = 1.4

mpl.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})


# List of HDF5 file names (e.g., "output_0001.h5", "output_0002.h5", etc.)
hdf5_filenames = [ "../hdf5/DWp_000036.2d.hdf5", "../hdf5/DWp_000061.2d.hdf5", "../hdf5/DWp_000068.2d.hdf5", "../hdf5/DWp_000100.2d.hdf5"]

def get_center(ds):
    # Small function that gets center for both single
    # or multiple dataset
    #   input: ds .. YT dataset
    #   return: center ... vector with center of box
    center = None
    if hasattr(ds,"domain_right_edge"):
        center = ds.domain_right_edge / 2.0
    elif hasattr(ds[0],"domain_right_edge"):
        center = ds[0].domain_right_edge / 2.0
    return center
center = get_center(yt.load(hdf5_filenames[0]))

width = [60.0,30.0]
center = [30.0, 15.0, 0]
center_array = []
for i, fn in enumerate(hdf5_filenames):
    center_array.append(center)
# Field to plot
def _field1(field, data):
        return data["rho"]  #* a_p**2

def plot_NG(a_data, ax, number):
    # number = str(a_data)[10:16]
    path = "../NG_a80_b40_prolate/a80_b40_prolate" + number + ".txt"

    dat = np.loadtxt(path) 
    x_values = abs(dat[:, 1]) - 30
    y_values = abs(dat[:, 0]) - 15


    ax.plot(x_values,y_values, linestyle = '--' ,color='k', linewidth = 2.5)

def plot_eik(a_data, ax, number):
 
     path = "../Eik_a59_b29/Eik_a59_b29_" + number + ".txt"

     dat = np.loadtxt(path) 
     x_values = abs(dat[:, 2]) - 30
     y_values = abs(dat[:, 0]) - 15


     ax.plot(x_values,y_values, linestyle = ':' ,color='k', linewidth = 3.5)

# Create a figure with 1 row and 4 columns
figure_size = (12, 12)
fig = plt.figure(figsize=figure_size)
grid = AxesGrid(fig, 111,
                    nrows_ncols = (2, 2),
                    axes_pad = 0.3,
                    label_mode = "L",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="7%",
                    cbar_pad="5%")

# Loop through each dataset, load it, and make a slice plot
for i, filename in enumerate(hdf5_filenames):
    # Load the HDF5 dataset
    ds = yt.load(filename)
    ds.add_field(("chombo","field1"), _field1, units="", sampling_type="cell")

    # Create the slice plot along the z-axis for the specified field
    slc = yt.SlicePlot(ds, "z", "field1", center = center_array[0], width=width, fontsize= F2)

    # Optional: Logarithmic scale for better visualization
    slc.set_log("field1", False)

    # Optional: Set colormap for the slice plot
    slc.set_cmap(field="field1", cmap="plasma")

    slc.set_zlim("field1", 1, 100.0)

   # slc.set_ylabel(r"$y$")
   # slc.set_xlabel(r"$x$")
    

    # Optional: Set units for axes
    #slc.set_axes_unit("kpc")

    # Render the plot as a matplotlib figure
    # slc.annotate_title(f"Time = {ds.current_time.in_units('Myr'):.1f} Myr")
   # fr = slc.to_mpl_figure()

    # Transfer the plot to the matplotlib axes
    plot = slc.plots["field1"]
    plot.figure = fig
    fig.figure_size = (12, 12)
    plot.axes = grid[i].axes
    grid[i].axes.set_xticks([])
    grid[i].axes.set_yticks([])
 #   grid[2].axes.set_xlabel(r"$y$")
   
    plot.cax = grid.cbar_axes[i]
    colorbar = plot.cb
    colorbar.fontsize= F2
    slc._setup_plots()

    xticks = [-30 ,-20, -10, 0, 10, 20, 30]
    x_labels = ['$0$', '$10$','$20$', '$30$','$40$', '$50$', '$60$']

    yticks = [-15, -5, 5, 15]
    y_labels = ['$0$', '$10$', '$20$', '$30$']

    
   

    # yticklabels = xticklabels

    grid[2].axes.set_xlabel(r"$x$",  fontsize=F2)
    grid[0].axes.set_xlabel("")
    grid[1].axes.set_xlabel("")

    grid[0].axes.set_ylabel(r"$y$",  fontsize=F2)
    grid[2].axes.set_ylabel(r"$y$",  fontsize=F2)
    grid[1].axes.set_xlabel("")
    grid[3].axes.set_xlabel("")
    time = ds.current_time
    grid[i].axes.set_title("$t = \:$%.1f" % time, fontsize=F2)
    plot.axes.set_xticks(xticks)
    plot.axes.set_yticks(yticks)
    plot.axes.set_xticklabels(x_labels)
    plot.axes.set_yticklabels(y_labels)

    plot.axes.tick_params(axis="both", direction="in",length=4, width=1,color="k",labelsize=F2)
    plot.axes.minorticks_off()
    plot.cax.set_ylabel("$\\rho$", rotation=0, fontsize=F2)
    time = ds.current_time
    grid[3].axes.set_xlabel(r"$x$",  fontsize=F2)
    plot_NG(ds, plot.axes, "%.1f" % time)
    t0 = abs(time.value - 27.2)
    if (time.value > 27.1):
        plot_eik(ds, plot.axes, "%.1f" % t0 )
    
  #  axes[i].set_title(f"Time = {ds.current_time.in_units('Myr'):.1f} Myr")
   # axes[i].axis('off')  # Hide axes for better visualization

# Adjust layout for tight spacing



plt.savefig("../plots/prolate_rho_plots_linear.png", dpi=200, bbox_inches = "tight")
plt.close()
