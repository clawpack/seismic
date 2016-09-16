

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

    def surface(current_data):
        """
        Return a masked array containing the surface elevation only in wet cells.
        Surface is q[2,:]
        """
        from numpy import ma
        drytol_default = rundata.geo_data.dry_tolerance
        drytol = getattr(current_data.user, 'dry_tolerance', drytol_default)
        q = current_data.q
        aux = current_data.aux
        h = q[0,:]
        eta = q[2,:]
        water = ma.masked_where(h<=drytol, eta)
        #water = eta
        return water

        
    def topo(current_data):
        q = current_data.q
        topo = q[2,:] - q[0,:]  # eta - h
        return topo

    def dtopo(current_data):
        q = current_data.q
        aux = current_data.aux
        topo = q[2,:] - q[0,:]  # eta - h
        topo_original = aux[2,:]
        dtopo = topo - topo_original
        return dtopo
        
    plotfigure = plotdata.new_plotfigure(name='domain', figno=0)
    plotfigure.kwargs = {'figsize': (8,9)}
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.xlimits = [-150e3,2e3]
    plotaxes.ylimits = [-.5, .5]
    plotaxes.title = 'Surface displacement'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = surface
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [-150e3,2e3]
    plotaxes.ylimits = [-3500, 500]
    plotaxes.title = 'Full depth'
    plotaxes.afteraxes = fixticks
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = surface
    plotitem.plot_var2 = topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [-150e3,2e3]
    plotaxes.ylimits = [-0.5,0.5]
    plotaxes.title = 'dtopo'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = dtopo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    #----------

    plotfigure = plotdata.new_plotfigure(name='shore', figno=1)
    #plotfigure.kwargs = {'figsize':(9,11)}


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [-50e3, 2000] #[180000,225200]
    plotaxes.ylimits = [-0.5,0.5]
    plotaxes.title = 'Zoom on shelf'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    #plotitem.plot_var = surface
    #plotitem.plot_var2 = topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = topo
    plotitem.color = 'k'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [-100,100]
    plotaxes.ylimits = [-2,2]
    plotaxes.title = 'Zoom around shore'

    plotaxes.afteraxes = fixticks

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False
    plotitem.plot_var = surface

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = surface
    plotitem.plot_var2 = topo
    plotitem.color = 'b'
    plotitem.MappedGrid = True
    plotitem.mapc2p = mapc2p

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = topo
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

