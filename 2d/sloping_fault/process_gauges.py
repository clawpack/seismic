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

subplot(212)
#ylim(-0.008,.008)
legend(loc='upper right', fontsize=8)

if plot_okada:
    from plot_okada import plot_okada_surface
    plot_okada_surface(ax, 'r-')
    legend()

