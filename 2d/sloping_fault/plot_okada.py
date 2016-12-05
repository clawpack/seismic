from pylab import *
from clawpack.clawutil.data import ClawData
import dtopotools_horiz_okada as dtopotools
reload(dtopotools)

# Read in fault info
column_map = {'mu':0,'dip':1,'width':2,'depth':3,'slip':4,'rake':5,'strike':6,
                'length':7,'longitude':8,'latitude':9,'rupture_time':10,'rise_time':11}
fault = dtopotools.Fault(coordinate_specification='top center')
fault.read('fault.data',column_map=column_map,skiprows=4)

for subfault in fault.subfaults:
    subfault.calculate_geometry()

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
