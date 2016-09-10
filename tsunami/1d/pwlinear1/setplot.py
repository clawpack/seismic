
import clawpack.geoclaw.shallow_1d.plot as geoplot

import setrun
rundata=setrun.setrun()

import mapc2p
reload(mapc2p)  # in case num_cells changed
from mapc2p import mapc2p

import numpy
from pylab import find


# Read in the maximum depth over the full run if available
# This will be plotted as a red line in the plots

try:
    fname = '_output/fort.hmax'
    d = numpy.loadtxt(fname)
    etamax = numpy.where(d[:,1]>1e-6, d[:,2], numpy.nan)
    xmax = d[:,0]
    jmax = find(d[:,1]>0).max()
    print "run-in = %8.2f,  run-up = %8.2f" % (d[jmax,0],d[jmax,2])
    print 'Loaded hmax from ',fname
except:
    xmax = None
    print "Failed to load fort.hmax"


def setplot(plotdata):

    plotdata.clearfigures()

    def fixticks1(current_data):
        from pylab import ticklabel_format
        ticklabel_format(format='plain',useOffset=False)

    def fixticks(current_data):
        from pylab import ticklabel_format, plot
        ticklabel_format(format='plain',useOffset=False)
        if xmax is not None:
            plot(xmax, etamax, 'r')

    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-200e3,2e3]
    plotaxes.ylimits = [-5, 5]
    plotaxes.title = 'Surface displacement'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [-200e3,2e3]
    plotaxes.ylimits = [-3500, 1000]
    plotaxes.title = 'Full depth'
    plotaxes.afteraxes = fixticks
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    #----------

    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    #plotfigure.kwargs = {'figsize':(9,11)}


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-50e3, 2000] #[180000,225200]
    plotaxes.ylimits = [-20, 20]
    plotaxes.title = 'Zoom on shelf'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.plot_var = geoplot.surface
    #plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    #plotaxes.xlimits = [-2000,2000]
    plotaxes.xlimits = [-1000,1000]
    #plotaxes.ylimits = [-10,40]
    plotaxes.ylimits = [-20,60]
    plotaxes.title = 'Zoom around shore'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = geoplot.surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = geoplot.surface
    plotitem.plot_var2 = geoplot.topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = geoplot.topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p


    plotdata.printfigs = True          # Whether to output figures
    plotdata.print_format = 'png'      # What type of output format
    plotdata.print_framenos = 'all'      # Which frames to output
    plotdata.print_fignos = 'all'      # Which figures to print
    plotdata.html = True               # Whether to create HTML files
    plotdata.latex = False             # Whether to make LaTeX output
    plotdata.parallel = True

    return plotdata

