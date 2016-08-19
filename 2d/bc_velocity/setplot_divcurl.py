
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import numpy as np

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps
    from matplotlib.colors import LinearSegmentedColormap 

    RGB = np.ones((100,3))
    for mred in range(10):
        for mblue in range(10):
            ind = mred*10 + mblue
            RGB[ind,0] = 1 - (mblue-1)/9.
            RGB[ind,1] = max(1 - (mred-1)/9.- (mblue-1)/9., 0)
            RGB[ind,2] = 1 - (mred-1)/9. 

    x = np.linspace(0,1,100)
    d = {}
    d['red'] = [(x[i],RGB[i,0],RGB[i,0]) for i in range(len(x))]
    d['green'] = [(x[i],RGB[i,1],RGB[i,1]) for i in range(len(x))]
    d['blue'] = [(x[i],RGB[i,2],RGB[i,2]) for i in range(len(x))]

    RXBmap = LinearSegmentedColormap('RXBmap',d)
    

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    def plot_interfaces(current_data):
        from pylab import linspace, plot
        xl = linspace(0,2,101)
        yl = 0.6 + 0.1*xl
        plot(xl,yl,'k')
        yl = 0.4 + 0.3*(xl-0.5) - 0.2*(xl-0.5)**2
        plot(xl,yl,'k')
    

    # Figure for sigma
    plotfigure = plotdata.new_plotfigure(name='trace(sigma)', figno=0)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'trace(sigma)'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces
    

    # Set up for item on these axes:

    def sigmatr(current_data):
        q = current_data.q
        return q[0,:,:] + q[1,:,:]

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = sigmatr
    plotitem.pcolor_cmap = colormaps.red_white_blue  
    plotitem.pcolor_cmin = -.01
    plotitem.pcolor_cmax = .01
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [False]
    

    def div(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            d = zeros(u.shape)
            return d
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        ux = (u[I+1,:][:,J] - u[I-1,:][:,J]) / (2*dx)
        vy = (v[:,J+1][I,:] - v[:,J-1][I,:]) / (2*dy)
        dint = ux + vy
        
        #zx = zeros((mx-2,1))
        #zy = zeros((1,my))
        #d = vstack((zy, hstack((zx, ux+vy, zx)), zy))
        
        d0 = dint[:,0]
        d1 = dint[:,-1]
        d2 = vstack((d0, dint.T, d1)).T
        d0 = d2[0,:]
        d1 = d2[-1,:]
        d = vstack((d0,d2,d1))      
        return d

    def curl(current_data):
        from numpy import array,zeros,hstack,vstack
        q = current_data.q
        u = q[3,:,:]
        v = q[4,:,:]
        mx, my = u.shape
        if (mx<3) or (my<3):
            c = zeros(u.shape)
            return c
        dx, dy = current_data.dx, current_data.dy
        I = array(range(1,mx-1))
        J = array(range(1,my-1))
        vx = (v[I+1,:][:,J] - v[I-1,:][:,J]) / (2*dx)
        uy = (u[:,J+1][I,:] - u[:,J-1][I,:]) / (2*dy)
        cint = vx - uy

        c0 = cint[:,0]
        c1 = cint[:,-1]
        c2 = vstack((c0, cint.T, c1)).T
        c0 = c2[0,:]
        c1 = c2[-1,:]
        c = vstack((c0,c2,c1))      
        return c
        
    def divcurl(current_data):
        from pylab import where, ceil
        dmax = 0.5
        cmax = 0.5
        d = abs(div(current_data)) / dmax
        c = abs(curl(current_data)) / cmax
        
        d = where(d<1, d, 1)
        c = where(c<1, c, 1)
        d = ceil(d*10)
        c = ceil(c*10)
        dc = d*10 + c
        #import pdb; pdb.set_trace()
        
        return dc
        
            
    # Figure for div and curl:
    # NOT WORKING!
    plotfigure = plotdata.new_plotfigure(name='PS', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'P- and S-waves'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = divcurl
    plotitem.pcolor_cmap = RXBmap
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 100
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [False]

    # Figure for div and curl side by side
    plotfigure = plotdata.new_plotfigure(name='div curl', figno=11)
    plotfigure.kwargs = {'figsize':(12,5)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(121)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'div(u) shows P-waves'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = div
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.3
    plotitem.pcolor_cmax = 0.3
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [1]

    #plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.plot_var = div
    #plotitem.contour_levels = np.linspace(-0.3,0.3,4)

    # Figure for curl:
    # plotfigure = plotdata.new_plotfigure(name='curl', figno=12)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(122)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'curl(u) shows S-waves'
    plotaxes.scaled = True
    plotaxes.afteraxes = plot_interfaces

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = curl
    plotitem.pcolor_cmap = colormaps.blue_white_red
    plotitem.pcolor_cmin = -0.3
    plotitem.pcolor_cmax = 0.3
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [False]
    plotitem.amr_patchedges_show = [1]

    if 0:
        # Set up for item on these axes:
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = False
        plotitem.plot_var = curl
        plotitem.contour_nlevels = 2
        plotitem.contour_min = 0.01
        plotitem.contour_max = .5
        plotitem.contour_colors = 'r'
        plotitem.amr_celledges_show = [False]
        plotitem.amr_patchedges_show = [False]



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Vertical velocity'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.plotstyle = 'b-'



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

    return plotdata

    
