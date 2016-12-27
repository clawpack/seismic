
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from mapping import Mapping
from clawpack.clawutil.data import ClawData
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    fault = dtopotools.Fault(coordinate_specification='top center')
    fault.read(plotdata.outdir + '/fault.data')

    mapping = Mapping(fault)
    fault_width = mapping.fault_width
    xcenter = mapping.xcenter
    ycenter = mapping.ycenter
    xp1 = mapping.xp1
    xp2 = mapping.xp2
    yp1 = mapping.yp1
    yp2 = mapping.yp2

    probdata = ClawData()
    probdata.read('setprob.data',force=True)
    xlimits = [xcenter-0.5*probdata.domain_width,xcenter+0.5*probdata.domain_width]
    ylimits = [-probdata.domain_depth,0.0]

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    def plot_fault(current_data):
        from pylab import linspace, plot
        xl = linspace(xp1,xp2,100)
        yl = linspace(yp1,yp2,100)
        plot(xl,yl,'k')

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:])

    def slip_direction_vel(current_data):
        # return vel dot tau, where tau is tangent to fault
        tau_x = (xp2 - xp1)/fault_width
        tau_y = (yp2 - yp1)/fault_width
        tau_x = np.where(current_data.y > ycenter, -tau_x, tau_x)
        tau_y = np.where(current_data.y > ycenter, -tau_y, tau_y)
        u = current_data.q[3,:,:]
        v = current_data.q[4,:,:]
        return u*tau_x + v*tau_y

    # Figure for waves
    plotfigure = plotdata.new_plotfigure(name='waves', figno=1)
    plotfigure.kwargs = {'figsize':(10,8)}

    # Set up axes for trace(sigma):
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
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
    plotitem.mapc2p = mapping.mapc2p

    # Set up axes for slip_direction_velocity:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'slip-direction-velocity'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_fault

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = slip_direction_vel
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = 'Level 6 grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#effeee', '#eeffee', '#eeeffe',
                                  '#eeeeff', '#ffffff']
    plotitem.amr_celledges_show = [0]
    plotitem.amr_patchedges_show = [0,0,0,0,0,1]
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapping.mapc2p

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
