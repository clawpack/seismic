from clawpack.clawutil.data import ClawData
import numpy
from pylab import *
from math import atan2

probdata = ClawData()
probdata.read('setprob.data',force=True)

fault_width = probdata.fault_width
theta = probdata.fault_dip
xcenter = probdata.fault_center
ycenter = -probdata.fault_depth
water_scaling = probdata.water_scaling

xcl = xcenter - 0.5*fault_width
xcr = xcenter + 0.5*fault_width

xp1 = xcenter - 0.5*fault_width*cos(theta)
xp2 = xcenter + 0.5*fault_width*cos(theta)
yp1 = ycenter - 0.5*fault_width*sin(theta)
yp2 = ycenter + 0.5*fault_width*sin(theta)

tol = min(abs(yp1),abs(yp2))

def mapc2p(xc,yc):
    """
    map computational grid to physical grid that rotates near the fault
    so cell edges match the fault line.  Linear interpolation is used to
    adjust the rotation angle based on distance from fault in computational space.
    The variable tol ensures the physical grid also lines up with a horizontal sea floor
    """

    # constucted signed distance function in computational domain
    ls = numpy.abs(yc - ycenter)
    ls = numpy.where(xc < xcl, numpy.sqrt((xc-xcl)**2 + (yc-ycenter)**2), ls)
    ls = numpy.where(xc > xcr, numpy.sqrt((xc-xcr)**2 + (yc-ycenter)**2), ls)

    # define grid that is rotated to line up with fault
    xrot = xcenter + numpy.cos(theta)*(xc-xcenter) - numpy.sin(theta)*(yc-ycenter)
    yrot = ycenter + numpy.sin(theta)*(xc-xcenter) + numpy.cos(theta)*(yc-ycenter)

    # Interpolate between roated grid and cartesian grid near the fault,
    # using cartesian grid far away from fault.
    xp = xc
    yp = numpy.where(yc > 0, yc*water_scaling, yc)
    xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xc, xp)
    yp = numpy.where(ls < tol, (tol-ls)/tol*yrot + ls/tol*yc, yp)

    return xp,yp

def test(mfault):

    domain_depth = probdata.domain_depth
    domain_width = probdata.domain_width
    fault_depth = probdata.fault_depth
    water_depth = probdata.water_depth

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    dx = fault_width/mfault
    mfault_to_floor = numpy.rint(fault_depth/dx)
    dy = fault_depth/mfault_to_floor
    mbelow_floor = int(numpy.ceil(domain_depth/dy))
    mabove_floor = int(numpy.floor(water_depth/dy))
    mx = int(numpy.ceil(domain_width/dx)) # mx
    my = mbelow_floor + mabove_floor# my
    mr = mx - mfault

    x = linspace(xcenter-0.5*fault_width - numpy.floor(mr/2.0)*dx, xcenter+0.5*fault_width + numpy.ceil(mr/2.0)*dx, mx+1)
    y = linspace(-mbelow_floor*dy, mabove_floor*dy, my+1)
    xc,yc = meshgrid(x,y)
    xp,yp = mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((xp1,xp2),(yp1,yp2),'-g', linewidth=2.0)
    plot(x,0*x,'-b', linewidth=2.0)

    axis('scaled')
