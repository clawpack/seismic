
from pylab import *

ylower = -50e3
yupper = 0.
yf1 = -15e3
yf2 = -23.6e3
#ycf = 1 + (yf1+yf2)/(2*(yupper-ylower))
ycf = 0.7

xlower = -50.e3
xupper = 150e3
xf1 = 0e3
xf2 = xf1 + 50e3

def mapc2p(xc,yc):
    """
    map computational rectangle (xlower,xupper) x (0,1) to 
    (xlower,xupper) x (ylower,yupper) with piecewise linear fault.
    The line yc=ycf maps to yp=yf1 where x<xf1 and to yp=yf2 where x>xf2,
    linear in between.
    """
    import numpy

    xp = xc
    yf = yf1 * numpy.ones(yc.shape)  # ok only where xc <= xf1
    yf = numpy.where(xc>xf1, yf1+(xc-xf1)*(yf2-yf1)/(xf2-xf1), yf)
    yf = numpy.where(xc>xf2, yf2, yf)

    yp = ylower + yc*(yf-ylower)/ycf
    yp = numpy.where(yc>ycf, yf + (yc-ycf)*(yupper-yf)/(1.-ycf), yp)

    return xp,yp
        
def test(mx,my):
    x = linspace(xlower,xupper,mx)
    y = linspace(0,1,my)
    xc,yc = meshgrid(x,y)
    xp,yp = mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
