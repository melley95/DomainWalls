import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
# import cmasher as cmr
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter
mpl.use("Agg")
import os
from unyt import unyt_array

F2 = 10 # Legend font size
line = 0.8 # Line width
alp = 1 # alpha
ticklabelsize = F2 

yt.enable_parallelism()

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

# plt.style.use("dark_background")


co = "w"
# ------------------- row 1 ---------------------

folder = "../hdf5/"
# folder = "/cosma7/data/dp092/dc-dejo1/k2/A0.0875_R1.4_Aalpha_0.5_w12/hdf5/"




def plot_NG(a_data, ax, number):
    # number = str(a_data)[10:16]
    path = "../NG_a80_b40_prolate/a80_b40_prolate" + number + ".txt"

    dat = np.loadtxt(path) 
    x_values = abs(dat[:, 1]) - 45
    y_values = abs(dat[:, 0]) - 25


    ax.plot(x_values,y_values, linestyle = '--' ,color='k', linewidth = 2.0)

def plot_eik(a_data, ax, number):
 
     path = "../Eik_a59_b29/Eik_a59_b29_" + number + ".txt"

     dat = np.loadtxt(path) 
     x_values = abs(dat[:, 2]) - 45
     y_values = abs(dat[:, 0]) - 25


     ax.plot(x_values,y_values, linestyle = ':' ,color='k', linewidth = 2.5)




# folder = "data/"
fns = [folder+"DWp_000000.2d.hdf5"
        ]

data_loc = folder + "DWp_*"
k=0

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

center = get_center(yt.load(fns[0]))

width = [90.0,50.0]
center = [45.0, 25.0, 0]
center_array = []
for i, fn in enumerate(fns):
    center_array.append(center)





ds_total = yt.load(data_loc)
for ds in ds_total[::1].piter():
    time = ds.current_time


  #  if(time > 70):
   #     L = 20
    #    width = [L,L,L, L,L,L, L,L,L]
    #    center = [10.0, 10.0, 0]
    #    center_array = []
    #    for i, fn in enumerate(fns):
    #        center_array.append(center)
# for ds in ds_total[::1]:
    # Load the data and create a single plot
    # ds = yt.load(fn) # load data
    # ds.length_unit.in_units("")
    all_data = ds.r[:,:]

 #  point = ds.r[0,0]
  #  a_p = point["scale_factor"].d
  #  print("scale", a_p)

    def _field1(field, data):
        return data["rho"]  #* a_p**2

 

    ds.add_field(("chombo","field1"), _field1, units="", sampling_type="cell")


    figure_size = (12, 8)
    fig = plt.figure(figsize=figure_size)
  


    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    # These choices of keyword arguments produce a four panel plot with a single
    # shared narrow colorbar on the right hand side of the multipanel plot. Axes
    # labels are drawn for all plots since we"re slicing along different directions
    # for each plot.
    grid = AxesGrid(fig, 121,
                    nrows_ncols = (1, 1),
                    axes_pad = 0.,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="5%")


  


    # fig.subplots_adjust(wspace=0.2, hspace=0)

    


    variable1 = "field1"

    p = yt.SlicePlot(ds, "z", variable1, center = center_array[0], width=width, fontsize= F2)
    
   
  

    p.set_log(variable1, True)


    # Ensure the colorbar limits match for all plots
    p.set_buff_size(2048)
   


    p.set_cmap(field=variable1, cmap="plasma")
    p.set_ylabel(r"$y$")
    p.set_xlabel(r"$x$")
    # p.hide_axes(draw_frame=True)

    max_var1 = all_data.max(variable1).d
    min_var1 = 0.0
    p.set_zlim(variable1, 0.0001, 10.0)




    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[variable1]
    plot.figure = fig
    fig.figure_size = (12, 8)
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]
    colorbar = plot.cb
    colorbar.fontsize= F2

    
    
    


    name = str(ds)[-12:-8]
    # print(ad["fieldphi"])
    # plot_AH(ds, p, name)
    # plot_AH(ds, p2, name)

    # Finally, this actually redraws the plot.
    p._setup_plots()




  #  x = np.linspace(-L/2., L/2., 1000)
    # plot.axes.axhline(y=0,xmin=0.5,lw=line,color=line_colors[i])
    # plot2.axes.axhline(y=0,xmin=0.5,lw=line,color=line_colors2[i])

    xticks = [-45 ,-25, -5, 15, 35]
    x_labels = ['$0$', '$20$', '$40$', '$60$', '$80$']

    yticks = [-25, -15, -5, 5, 15, 25]
    y_labels = ['$0$', '$10$', '$20$', '$30$', '$40$', '$50$']

   

   

    # yticklabels = xticklabels

    plot.axes.set_xticks(xticks)
    plot.axes.set_yticks(yticks)
    plot.axes.set_xticklabels(x_labels)
    plot.axes.set_yticklabels(y_labels)
  
  #  plot.axes.minorticks_off()
    plot.axes.tick_params(axis="both", direction="in",length=4, width=1,color=co,labelsize=ticklabelsize)
 #  plot.title("\phi")
   
    plot.axes.set_title(r"$t=%.1f$"%time,c="k",fontsize=F2)
    plot.cax.set_ylabel("$\\rho$", rotation=0, fontsize=F2)




    from matplotlib.ticker import ScalarFormatter
    class ScalarFormatterClass(ScalarFormatter):
        def _set_format(self):
            self.format = "$%1.1f$"
    yScalarFormatter = ScalarFormatterClass(useMathText=True)
    yScalarFormatter.set_powerlimits((0,0))


    dvar1 = 0.2*(max_var1 - min_var1)
    plot.axes.spines['top'].set_visible(False)
    plot.axes.spines['right'].set_visible(False)




    plot.axes.tick_params(top=False, right=False)

    plot.axes.minorticks_off()


   # from matplotlib.ticker import ScalarFormatter
   # class ScalarFormatterClass(ScalarFormatter):
    #    def _set_format(self):
     #       self.format = "$%1.1f$"
  #  yScalarFormatter = ScalarFormatterClass(useMathText=True)
  #  yScalarFormatter.set_powerlimits((0,0))

    mintick1 = min_var1+dvar1
    maxtick1 = max_var1-dvar1





    plot_NG(ds, plot.axes, "%.1f" % time)
    t0 = abs(time.value - 27.2)
    if (time.value > 27.1):
        plot_eik(ds, plot.axes, "%.1f" % t0 )
    
    
    plot.axes.legend(['NG','Eikonal'],fontsize=0.8*F2)
 
    # fig.suptitle(r"$t=%i$"%int(time),c="k",fontsize=F2*1.2, x=0.3, y=0.7, ha='center', va='top')

    # plt.savefig("frames_movie/"+"movie_z%04i.png"%k, dpi=200, bbox_inches = "tight")
    plt.savefig("rho_movie/"+"movie_z"+name+".png", dpi=200, bbox_inches = "tight")
    k += 1
    # print(i)
    plt.close()