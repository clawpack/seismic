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

figure(500)
clf()

ngauges = 100

for t in linspace(0.,60,7):
    xs = zeros(ngauges)
    ys = zeros(ngauges)
    for j in range(ngauges):
        gaugeno = j            # at top of domain
        #gaugeno = 200 + j      # at bottom of domain
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

subplot(211)
#ylim(-0.008,.008)
legend(loc='lower right', fontsize=10)

