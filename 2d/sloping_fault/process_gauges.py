from pylab import *
from clawpack.visclaw.data import ClawPlotData

gdata = loadtxt('gauges.data',skiprows=7)
ngauges = gdata.shape[0]
print "Found %s gauges" % ngauges
xc = gdata[:,1]
yc = gdata[:,2]

# Load ClawplotData
plotdata = ClawPlotData()
plotdata.outdir = '_output'

figure(501)
clf()

plot_okada = False
ngauges = 100;  goffset = 0; plot_okada = True  # top surface
#ngauges = 100;  goffset = 200  # bottom
#ngauges = 50;  goffset = 300  # above fault plane
#ngauges = 50;  goffset = 400  # below fault plane

#for t in linspace(0.,120,7):
for t in [100]:   #linspace(0.,100,11):
    xs = zeros(ngauges)
    ys = zeros(ngauges)
    for j in range(ngauges):
        gaugeno = j            # at top of domain
        gaugeno = goffset + j      # at bottom of domain
        g = plotdata.getgauge(gaugeno)
        for k in range(1,len(g.t)):
            if g.t[k] > t:
                break
            dt = g.t[k] - g.t[k-1]
            u = g.q[3,k]
            v = g.q[4,k]
            xs[j] = xs[j] + dt*u
            ys[j] = ys[j] + dt*v

    subplot(211)
    plot(xc[:ngauges],xs,label='t = %6.3f' % t)
    title("horizontal displacement")
    subplot(212)
    plot(xc[:ngauges],ys,label='t = %6.3f' % t)
    title("vertical displacement")

subplot(212)
#ylim(-0.008,.008)
legend(loc='upper right', fontsize=8)

if plot_okada:
    from clawpack.geoclaw import dtopotools
    fault = dtopotools.Fault()
    subfault = dtopotools.SubFault()
    subfault.mu = 3e10
    subfault.depth = -15e3
    subfault.dip = 9.9
    subfault.width = sqrt(50e3**2 + 8.6e3**2)
    subfault.slip = 1.
    subfault.rake = 90
    subfault.strike = 0
    subfault.length = 1000e3
    subfault.longitude = 0.
    subfault.latitude = 0.
    subfault.coordinate_specification = 'top center'
    subfault.calculate_geometry()

    fault.subfaults = [subfault] 
    print "Mw = ", fault.Mw()

    x = linspace(-1.0,1.5,1000)
    y = array([0.])
    dtopo = fault.create_dtopography(x,y,[1.])

    plot(dtopo.x*111e3,dtopo.dZ[0,0,:],'k--')

