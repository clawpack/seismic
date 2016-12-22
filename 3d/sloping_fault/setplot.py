
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import os, shutil
from mapc2p import mapc2p
from clawpack.clawutil.data import ClawData
cscale = 8 # scale color limits

probdata = ClawData()
probdata.read('setprob.data',force=True)

width = probdata.fault_width
theta = probdata.fault_dip
xcenter = probdata.fault_xcenter
zcenter = -probdata.fault_depth

xp1 = xcenter - 0.5*width*np.cos(theta)
xp2 = xcenter + 0.5*width*np.cos(theta)
zp1 = zcenter + 0.5*width*np.sin(theta)
zp2 = zcenter - 0.5*width*np.sin(theta)

xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
zlimits = [-probdata.domain_depth,0.0]

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """
    slice_number = 3
    os.chdir(plotdata.outdir)
    for filename in os.listdir('.'):
        if (filename.startswith('slice_%d' % slice_number)):
            shutil.copyfile(filename,filename.replace('slice_%d' % slice_number,'fort',1))

    os.chdir('..')


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def plot_interfaces(current_data):
        from pylab import linspace, plot, sin, pi
        xl = linspace(xp1,xp2,100)
        zl = linspace(zp1,zp2,100)
        #yl = -8e3 - sin(10*pi/180.)*xl
        plot(xl,zl,'k')


    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])

    # Figure for trace(sigma)
    plotfigure = plotdata.new_plotfigure(name='trace', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = zlimits
    plotaxes.title = '-trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e6
    plotitem.pcolor_cmax = 1e6
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = zlimits
    plotaxes.title = 'z-velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 8
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = zlimits
    plotaxes.title = 'Level 4 grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
#    plotdata.parallel = True

    return plotdata
