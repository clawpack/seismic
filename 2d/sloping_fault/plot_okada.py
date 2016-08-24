from pylab import *
from clawpack.clawutil.data import ClawData
import dtopotools_horiz_okada as dtopotools
reload(dtopotools)
fault = dtopotools.Fault()
fault.subfaults = []

# Read in fault info
probdata = ClawData()
probdata.read('setprob.data',force=True)

width = probdata.fault_width
theta = probdata.fault_dip
xcenter = probdata.fault_center
ycenter = -probdata.fault_depth
domain_width = probdata.domain_width
ndip = 50

dlongitude = (width/111.e3) / ndip   # convert to degees and split up
sf_width = width/ndip
x = arange(xcenter-0.5*width,xcenter+0.5*width,sf_width)
slip = exp(-((x-xcenter)/(0.5*width))**2)

for i in range(ndip):
    subfault = dtopotools.SubFault()
    subfault.mu = 3e10
    subfault.dip = theta/pi*180.0
    subfault.width = sf_width
    subfault.depth = -(ycenter + (x[i]-xcenter)*sin(theta))
    subfault.slip = slip[i]
    subfault.rake = 90
    subfault.strike = 0
    subfault.length = 1000e3
    subfault.longitude = dlongitude*i
    subfault.latitude = 0.
    subfault.coordinate_specification = 'top center'
    subfault.calculate_geometry()

    fault.subfaults.append(subfault)

print "Mw = ", fault.Mw()

x = linspace(-1.5,1.5,1000)
y = array([0.])
dtopo = fault.create_dtopography(x,y,[1.],y_disp=True)

def plot_okada_surface(ax=None, plotstyle='k-'):
    if ax is None:
        figure()
        ax = subplot(111)
    ax.plot(dtopo.x*111.e3,dtopo.dZ[0,0,:],plotstyle,label="Okada")

def plot_okada_horiz(ax=None, plotstyle='k-'):
    if ax is None:
        figure()
        ax = subplot(111)
    ax.plot(dtopo.x*111.e3,dtopo.dY[0,0,:],plotstyle,label="Okada")

if __name__=='__main__':
    figure()
    ax = subplot(211)
    plot_okada_surface(ax)
    ax = subplot(212)
    plot_okada_horiz(ax)
    savefig('okada.png')
    print "Created okada.png"
