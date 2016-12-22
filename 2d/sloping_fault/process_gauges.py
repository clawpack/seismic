from pylab import *
from clawpack.visclaw.data import ClawPlotData


if 0:
    gdata = loadtxt('gauges.data',skiprows=7)
    ngauges = gdata.shape[0]
    print "Found %s gauges" % ngauges
    xc = gdata[:,1]
    yc = gdata[:,2]

# Load ClawplotData
plotdata = ClawPlotData()
plotdata.format = 'binary'
plotdata.outdir = '_output'



def plot_gauges(goffset=0):

    if goffset == 0:
        ngauges = 100 # top surface
    elif goffset == 200:
        ngauges = 100 # bottom of domain
    elif goffset == 300:
        ngauges = 50 # above fault plane
    elif goffset == 400:
        ngauges = 50 # below fault plane
    plot_okada = (goffset==0)
    
    figure()

    #for t in linspace(0.,120,7):
    for t in [100]:   #linspace(0.,100,11):
        xs = zeros(ngauges)
        ys = zeros(ngauges)
        xg = zeros(ngauges)
        for j in range(ngauges):
            gaugeno = goffset + j
            g = plotdata.getgauge(gaugeno)
            xg[j] = g.location[0]  # x-location of this gauge
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                u = g.q[3,k]
                v = g.q[4,k]
                xs[j] = xs[j] + dt*u
                ys[j] = ys[j] + dt*v
    
        subplot(211)
        plot(xg,xs,label='t = %6.3f' % t)
        title("horizontal displacement")
        ax = subplot(212)
        plot(xg,ys,label='t = %6.3f' % t)
        title("vertical displacement")
    
    if plot_okada:
        import plot_okada
        ax = subplot(211)
        plot_okada.plot_okada_horiz(ax, 'r-')
        legend()
    
        ax = subplot(212)
        plot_okada.plot_okada_surface(ax, 'r-')
        legend()


def fault_slip(t):
    dip = 9.9  # degrees
    ngauges = 50
    xs_above = zeros(ngauges)
    ys_above = zeros(ngauges)
    xs_below = zeros(ngauges)
    ys_below = zeros(ngauges)
    xg = zeros(ngauges)
    slip = zeros(ngauges)
    for j in range(ngauges):
        g_above = plotdata.getgauge(j + 300)
        g_below = plotdata.getgauge(j + 400)
        xg[j] = g_above.location[0]  # x-location of this gauge
        for k in range(1,len(g_above.t)):
            if g_above.t[k] > t:
                break
            dt = g_above.t[k] - g_above.t[k-1]
            u_above = g_above.q[3,k]
            v_above = g_above.q[4,k]
            u_below = g_below.q[3,k]
            v_below = g_below.q[4,k]
            xs_above[j] = xs_above[j] + dt*u_above
            ys_above[j] = ys_above[j] + dt*v_above
            xs_below[j] = xs_below[j] + dt*u_below
            ys_below[j] = ys_below[j] + dt*v_below
            
    dipdir_above = xs_above*cos(dip*pi/180.) - ys_above*sin(dip*pi/180.)
    orthog_above = xs_above*sin(dip*pi/180.) + ys_above*cos(dip*pi/180.)
    dipdir_below = xs_below*cos(dip*pi/180.) - ys_below*sin(dip*pi/180.)
    orthog_below = xs_below*sin(dip*pi/180.) + ys_below*cos(dip*pi/180.)
    slip = dipdir_below - dipdir_above
    return xg, slip


def make_dtopo(xtopo, times):
    """
    Parse the gauge output for gauges on surface (seafloor)
    and turn into dtopo object.  Interpolate to the requested xtopo, times.
    """

    from scipy.interpolate import interp1d

    # custom dtopotools includes DTopography1d
    import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools

    goffset = 0  # gauges at top surface
    ngauges = 100

    # dtopo object for vertical displacement:
    dtopo_z = dtopotools.DTopography1d()
    dtopo_z.x = xtopo
    dtopo_z.times = times

    # dtopo object for horizontal displacement:
    dtopo_x = dtopotools.DTopography1d()
    dtopo_x.x = xtopo
    dtopo_x.times = times

    dZ_list_z = []
    dZ_list_x = []

    for t in times:
        xs = zeros(ngauges)
        ys = zeros(ngauges)
        xg = zeros(ngauges)
        for j in range(ngauges):
            gaugeno = goffset + j
            g = plotdata.getgauge(gaugeno)
            xg[j] = g.location[0]  # x-location of this gauge
            for k in range(1,len(g.t)):
                if g.t[k] > t:
                    break
                dt = g.t[k] - g.t[k-1]
                u = g.q[3,k]
                v = g.q[4,k]
                xs[j] = xs[j] + dt*u
                ys[j] = ys[j] + dt*v

        # vertical displacement:
        dz = interp1d(xg, ys, 'linear',bounds_error=False,
                fill_value=0.,assume_sorted=True)
        dZ_list_z.append(dz(xtopo))

        # horizontal displacement:
        dz = interp1d(xg, xs, 'linear',bounds_error=False,
                fill_value=0.,assume_sorted=True)
        dZ_list_x.append(dz(xtopo))

    dtopo_z.dZ = array(dZ_list_z, ndmin=2)
    dtopo_x.dZ = array(dZ_list_x, ndmin=2)

    return dtopo_z, dtopo_x


def save_dtopo_test1():
    """
    Test make_dtopo.
    Save full time history and also make a dtopo file with the final 
    deformation shifted to t = 1 second.
    """
    
    import copy
    xtopo = linspace(-150e3,150e3,151)
    times = linspace(0,80,41)
    dtopo_z, dtopo_x = make_dtopo(xtopo,times)

    fname = 'dtopo_z_seismic.tt3'
    dtopo_z.write(fname, 3)
    print "Created ",fname

    fname = 'dtopo_x_seismic.tt3'
    dtopo_x.write(fname, 3)
    print "Created ",fname

    # make the deformation file with only the first and last columns of
    # dtopo.dZ (initial and final time) and set final time to 1 sec:
    dtopo_final = copy.copy(dtopo_z)
    dtopo_final.times = array([0., 1.])
    dtopo_final.dZ = dtopo_z.dZ[[0,-1],:]
    fname = 'dtopo_z_seismic_final.tt3'
    dtopo_final.write(fname, 3)
    print "Created ",fname

    dtopo_final = copy.copy(dtopo_x)
    dtopo_final.times = array([0., 1.])
    dtopo_final.dZ = dtopo_x.dZ[[0,-1],:]
    fname = 'dtopo_x_seismic_final.tt3'
    dtopo_final.write(fname, 3)
    print "Created ",fname

