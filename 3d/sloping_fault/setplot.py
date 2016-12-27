
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

    fault = dtopotools.Fault(coordinate_specification='top center')
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
    zlimits = [-probdata.domain_depth,0.0]

    def plot_fault_xz(current_data):
        from pylab import linspace, plot
        xl = linspace(xp1,xp2,100)
        zl = linspace(zp1,zp2,100)
        plot(xl,zl,'k')

    if (slice_number is 1):
        x1limits = xlimits
        x2limits = ylimits
        mapc2p = mapping.mapc2p_xy
        plot_fault = None
    elif (slice_number is 2):
        x1limits = ylimits
        x2limits = zlimits
        mapping.set_slice_xval(0.0)
        mapc2p = mapping.mapc2p_yz
        plot_fault = None
    elif (slice_number is 3):
        x1limits = xlimits
        x2limits = zlimits
        mapc2p = mapping.mapc2p_xz
        plot_fault = plot_fault_xz

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:] + q[2,:,:])

    # Figure for trace(sigma)
    plotfigure = plotdata.new_plotfigure(name='trace', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = x1limits
    plotaxes.ylimits = x2limits
    plotaxes.title = '-trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

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
    plotaxes.xlimits = x1limits
    plotaxes.ylimits = x2limits
    plotaxes.title = 'x-velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 6
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -1e-1
    plotitem.pcolor_cmax = 1e-1
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = x1limits
    plotaxes.ylimits = x2limits
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
