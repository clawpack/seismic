from pylab import *
from clawpack.clawutil.data import ClawData
import clawpack.geoclaw.dtopotools as dtopotools
#reload(dtopotools)
fault = dtopotools.Fault()
fault.subfaults = []

# Read in fault info
probdata = ClawData()
probdata.read('setprob.data',force=True)

width = probdata.fault_width
length = probdata.fault_length
theta = probdata.fault_dip
xcenter = probdata.fault_xcenter
ycenter = probdata.fault_ycenter
zcenter = probdata.fault_depth
domain_width = probdata.domain_width
domain_length = probdata.domain_length
ndip = 30
nstrike = 30

dlongitude = (width/111.e3) / ndip   # convert to degees and split up
dlatitude = (length/111.e3) / nstrike

sf_width = width/ndip
sf_length = length/ndip
x = arange(xcenter-0.5*width,xcenter+0.5*width,sf_width)
y = arange(ycenter-0.5*length,ycenter+0.5*length,sf_length)
X,Y = meshgrid(x,y)
slip = exp(-((X-xcenter)/(0.5*width))**2 - ((Y-ycenter)/(0.5*length))**2)

for i in range(ndip):
    for j in range(nstrike):
        subfault = dtopotools.SubFault()
        subfault.mu = 3e10
        subfault.dip = theta/pi*180.0
        subfault.width = sf_width
        subfault.depth = zcenter + (x[i]-xcenter)*sin(theta)
        subfault.slip = slip[i,j]
        subfault.rake = 90
        subfault.strike = 0
        subfault.length = sf_length
        subfault.longitude = dlongitude*i
        subfault.latitude = dlatitude*j
        subfault.coordinate_specification = 'top center'
        subfault.calculate_geometry()

        fault.subfaults.append(subfault)

print "Mw = ", fault.Mw()

x = linspace(-1.5,1.5,100)
y = linspace(-1.5,1.5,100)
z = array([0.])
dtopo = fault.create_dtopography(x,y,[1.])

def plot_okada(ax=None):
    if ax is None:
        figure()
        ax = subplot(111)
#    dtopotools.plot_dZ_colors(dtopo.x*111.e3,dtopo.y*111.e3,dtopo.dZ[0,:,:],axes=ax)
    pcolor(dtopo.x*111.e3,dtopo.y*111.e3,dtopo.dZ[0,:,:],axes=ax)

if __name__=='__main__':
    dtopotools.plot_dZ_colors(dtopo.x*111.e3,dtopo.y*111.e3,dtopo.dZ[0,:,:])
    savefig('okada.png')
    print "Created okada.png"
