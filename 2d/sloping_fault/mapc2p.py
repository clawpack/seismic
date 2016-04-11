from clawpack.clawutil.data import ClawData
import numpy
from pylab import *
from math import atan2

probdata = ClawData()
probdata.read('setprob.data',force=True)

width = probdata.fault_width
theta = probdata.fault_dip
xcenter = probdata.fault_center
ycenter = -probdata.fault_depth

xcl = xcenter - 0.5*width
xcr = xcenter + 0.5*width

xp1 = xcenter - 0.5*width*cos(theta)
xp2 = xcenter + 0.5*width*cos(theta)
yp1 = ycenter - 0.5*width*sin(theta)
yp2 = ycenter + 0.5*width*sin(theta)

tol = min(abs(yp1),abs(yp2))

def mapc2p(xc,yc):
    """
    map computational rectangle (xlower,xupper) x (0,1) to
    (xlower,xupper) x (ylower,yupper) with piecewise linear fault.
    The line yc=ycf maps to yp=yf1 where x<xf1 and to yp=yf2 where x>xf2,
    linear in between.
    """

    ls = numpy.abs(yc - ycenter)
    ls = numpy.where(xc < xcl, numpy.sqrt((xc-xcl)**2 + (yc-ycenter)**2), ls)
    ls = numpy.where(xc > xcr, numpy.sqrt((xc-xcr)**2 + (yc-ycenter)**2), ls)

    xrot = xcenter + numpy.cos(theta)*(xc-xcenter) - numpy.sin(theta)*(yc-ycenter)
    yrot = ycenter + numpy.sin(theta)*(xc-xcenter) + numpy.cos(theta)*(yc-ycenter)

    xp = xc
    yp = yc
    xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xc, xp)
    yp = numpy.where(ls < tol, (tol-ls)/tol*yrot + ls/tol*yc, yp)

    return xp,yp

def test(mx,my):
    ratio = probdata.fault_depth/probdata.domain_depth
    num_cells_above = numpy.ceil(ratio*my)
    y = linspace(0.0,probdata.fault_depth,num_cells_above)
    y = numpy.append(y,linspace(probdata.fault_depth,probdata.domain_depth,my - num_cells_above))
    y = -y

    x = linspace(-0.5*probdata.domain_width,0.5*probdata.domain_width,mx)
    xc,yc = meshgrid(x,y)
    xp,yp = mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((xp1,xp2),(yp1,yp2),'-g')
    axis('scaled')
