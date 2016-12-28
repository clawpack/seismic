
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import os, shutil
from mapping import Mapping
from clawpack.clawutil.data import ClawData
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from clawpack.geoclaw.data import LAT2METER

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """
    slice_number = 1 # set to surface slice number
    os.chdir(plotdata.outdir)
    for filename in os.listdir('.'):
        if (filename.startswith('slice_%d' % slice_number)):
            shutil.copyfile(filename,filename.replace('slice_%d' % slice_number,'fort',1))

    os.chdir('..')

    fault = dtopotools.Fault()
    fault.read('fault.data')

    mapping = Mapping(fault)

    xp1 = mapping.xp1
    xp2 = mapping.xp2
    zp1 = mapping.zp1
    zp2 = mapping.zp2
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter

    probdata = ClawData()
    probdata.read('setprob.data',force=True)
    probdata.read(plotdata.outdir + '/setprob.data',force=True)
    xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
    ylimits = [ycenter-0.5*probdata.domain_length,ycenter+0.5*probdata.domain_length]

    from clawpack.visclaw import colormaps

    xc = np.linspace(-150e3,200e3,350)
    yc = np.linspace(-87.5e3,87.5e3,175)
    fault.create_dtopography(xc/LAT2METER,yc/LAT2METER.e3,[1.0])

    okada_max = np.max(fault.dtopo.dZ[0,:,:])
    clevels = np.linspace(-okada_max,okada_max,10)

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    def plot_okada_contour(current_data):
        from pylab  import gca

        kwargs = {'levels':clevels,'colors':'g'}
        ax = gca()
        fault.plot_okada_contour(axes=ax,kwargs=kwargs)

    def plot_okada(current_data):
        from pylab import gca

        kwargs = {'cmap':colormaps.blue_white_red,'vmin':clevels[0],'vmax':clevels[-1]}
        ax = gca()
        fault.plot_okada(axes=ax,dim=2,kwargs=kwargs)
        kwargs = {'levels':clevels,'colors':'g','linewidths':3}
        fault.plot_okada_contour(axes=ax,kwargs=kwargs)

    # Figure for vertical displacement
    plotfigure = plotdata.new_plotfigure(name='trace', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set axes for numerical solution:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Numerical solution'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_okada_contour

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 9
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = clevels[0]
    plotitem.pcolor_cmax = clevels[-1]
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 9
    plotitem.contour_colors = 'k'
    plotitem.contour_levels = clevels
    plotitem.kwargs = {'linewidths':3}
    plotitem.amr_contour_show = [0,0,0,1]
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]

    # Set axes for Okada solution:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title_with_t = False
    plotaxes.title = 'Okada solution'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_okada

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Level 4 grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,1]

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
