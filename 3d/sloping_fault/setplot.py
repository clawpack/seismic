
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
from plot_okada import plot_okada
from clawpack.clawutil.data import ClawData


cscale = 8 # scale color limits

probdata = ClawData()
probdata.read('setprob.data',force=True)
xlimits = [-0.5*probdata.domain_width,0.5*probdata.domain_width]
ylimits = [-0.5*probdata.domain_length,0.5*probdata.domain_length]

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def afterframe(current_data):
        from pylab import figure,subplot,linspace,title,zeros,pcolor,colorbar,axis
        from clawpack.visclaw.data import ClawPlotData
        ngauges = 100;  goffset = 0
        t = current_data.t

        x = zeros((ngauges,ngauges))
        y = zeros((ngauges,ngauges))
        dz = zeros((ngauges,ngauges))
        for i in range(ngauges):
            for j in range(ngauges):
                gaugeno = goffset + i*ngauges + j
                g = plotdata.getgauge(gaugeno,verbose=False)
                for k in range(1,len(g.t)):
                    if g.t[k] > t:
                        break
                    dt = g.t[k] - g.t[k-1]
                    w = g.q[8,k]
                    x[i,j] = g.location[0]
                    y[i,j] = g.location[1]
                    dz[i,j] = dz[i,j] + dt*w

        figure(10)
        ax = subplot(3,1,1)
        pcolor(x,y,dz,cmap='RdBu')
        axis([xlimits[0], xlimits[1], ylimits[0], ylimits[1]])
    #    colorbar()
        title("surface displacement")
        ax = subplot(3,1,2)
        plot_okada(ax)
        axis([xlimits[0], xlimits[1], ylimits[0], ylimits[1]])

    plotdata.afterframe = afterframe

    def sigmatr(current_data):
        # return -trace(sigma)
        q = current_data.q
        return -(q[0,:,:] + q[1,:,:] + q[2,:,:])

    # Figure for trace(sigma)
    plotfigure = plotdata.new_plotfigure(name='trace', figno=10)
    plotfigure.kwargs = {'figsize':(8,8)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,3)'
    #plotaxes.xlimits = [-75e3, 125e3]
    #plotaxes.ylimits = [-50e3,0]
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = '-trace(sigma)'
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -2e6
    plotitem.pcolor_cmax = 2e6
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0]
    plotitem.amr_patchedges_show = [0]
    plotitem.MappedGrid = False
#    plotitem.mapc2p = mapc2p




    # # Figure for trace(sigma) and sigma_12 side by side
    # plotfigure = plotdata.new_plotfigure(name='P and S waves', figno=11)
    # plotfigure.show = False
    # plotfigure.kwargs = {'figsize':(12,12)}
    #
    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(511)'
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = '-trace(sigma)'
    # plotaxes.scaled = True
    # #plotaxes.afteraxes = plot_interfaces
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = sigmatr
    # plotitem.pcolor_cmap = colormaps.blue_white_red
    # plotitem.pcolor_cmin = -0.03 * cscale
    # plotitem.pcolor_cmax = 0.03 * cscale
    # plotitem.add_colorbar = False
    # plotitem.amr_celledges_show = [1,0]
    # plotitem.amr_patchedges_show = [0]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p
    #
    #
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(512)'
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'div(u)'
    # plotaxes.scaled = True
    # #plotaxes.afteraxes = plot_interfaces
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = div
    # plotitem.pcolor_cmap = colormaps.blue_white_red
    # plotitem.pcolor_cmin = -0.1 * cscale
    # plotitem.pcolor_cmax = 0.1 * cscale
    # plotitem.add_colorbar = False
    # plotitem.amr_celledges_show = [False]
    # plotitem.amr_patchedges_show = [0]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p
    #
    #
    # # Figure for curl:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(513)'
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'curl(u)'
    # plotaxes.scaled = True
    # #plotaxes.afteraxes = plot_interfaces
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = curl
    # plotitem.pcolor_cmap = colormaps.blue_white_red
    # plotitem.pcolor_cmin = -0.2 * cscale
    # plotitem.pcolor_cmax = 0.2 * cscale
    # plotitem.add_colorbar = False
    # plotitem.colorbar_shrink = 0.7
    # plotitem.amr_celledges_show = [False]
    # plotitem.amr_patchedges_show = [0]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p
    #
    #
    # # Figure for x-velocity:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(514)'
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'x-velocity'
    # plotaxes.scaled = True
    # #plotaxes.afteraxes = plot_interfaces
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = 3
    # plotitem.pcolor_cmap = colormaps.blue_white_red
    # plotitem.pcolor_cmin = -0.005 * cscale
    # plotitem.pcolor_cmax = 0.005 * cscale
    # plotitem.add_colorbar = False
    # plotitem.colorbar_shrink = 0.7
    # plotitem.amr_celledges_show = [False]
    # plotitem.amr_patchedges_show = [0]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p
    #
    #
    # # Figure for y-velocity:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(515)'
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'y-velocity'
    # plotaxes.scaled = True
    # #plotaxes.afteraxes = plot_interfaces
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = 4
    # plotitem.pcolor_cmap = colormaps.blue_white_red
    # plotitem.pcolor_cmin = -0.005 * cscale
    # plotitem.pcolor_cmax = 0.005 * cscale
    # plotitem.add_colorbar = False
    # plotitem.colorbar_shrink = 0.7
    # plotitem.amr_celledges_show = [False]
    # plotitem.amr_patchedges_show = [0]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p
    #


    # # Figure for grid cells
    # plotfigure = plotdata.new_plotfigure(name='cells', figno=2)
    #
    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = xlimits
    # plotaxes.ylimits = ylimits
    # plotaxes.title = 'Level 4 grid patches'
    # plotaxes.scaled = True
    #
    # # Set up for item on these axes:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    # plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee', '#ffffff']
    # plotitem.amr_celledges_show = [0]
    # plotitem.amr_patchedges_show = [0,0,0,1]
    # plotitem.MappedGrid = True
    # plotitem.mapc2p = mapc2p


    # #-----------------------------------------
    # # Figures for gauges
    # #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
    #                 type='each_gauge')
    # #plotfigure.clf_each_gauge = False
    # plotfigure.show = False
    #
    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(2,1,1)'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'Horizontal velocity'
    #
    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 3
    # plotitem.plotstyle = 'b-'
    #
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(2,1,2)'
    # plotaxes.ylimits = 'auto'
    # plotaxes.title = 'Vertical velocity'
    #
    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 4
    # plotitem.plotstyle = 'b-'



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
    plotdata.parallel = True

    return plotdata
