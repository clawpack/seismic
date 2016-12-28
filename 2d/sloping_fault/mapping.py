import numpy
from pylab import *
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)
from clawpack.geoclaw.data import LAT2METER

def test(mfault):

    from clawpack.clawutil.data import ClawData

    probdata = ClawData()
    probdata.read('setprob.data',force=True)

    fault = dtopotools.Fault()
    fault.read('fault.data')

    mapping = Mapping(fault)

    domain_depth = probdata.domain_depth
    domain_width = probdata.domain_width

    # num of cells here determined in an identical fashion to that in setrun.py
    # additional comments can be found there
    dx = mapping.fault_width/mfault
    num_cells_above = numpy.rint(mapping.fault_depth/dx)
    dy = mapping.fault_depth/num_cells_above
    mx = int(numpy.ceil(domain_width/dx)) # mx
    my = int(numpy.ceil(domain_depth/dy)) # my
    mr = mx - mfault

    x = linspace(mapping.xcenter-0.5*mapping.fault_width - numpy.floor(mr/2.0)*dx, mapping.xcenter+0.5*mapping.fault_width + numpy.ceil(mr/2.0)*dx, mx+1)
    y = linspace(-my*dy, 0.0, my+1)
    xc,yc = meshgrid(x,y)
    xp,yp = mapping.mapc2p(xc,yc)
    figure()
    plot(xp,yp,'k-')
    plot(xp.T,yp.T,'k-')
    plot((mapping.xp1,mapping.xp2),(mapping.yp1,mapping.yp2),'-g')
    axis('scaled')


class Mapping(object):

    def __init__(self, fault):

        subfaultF = fault.subfaults[0]
        subfaultL = fault.subfaults[-1]
        theta = subfaultL.dip/180.0*numpy.pi

        xp1 = subfaultF.longitude*LAT2METER
        yp1 = -subfaultF.depth
 
        xp2 = subfaultL.longitude*LAT2METER + np.cos(theta)*subfaultL.width
        yp2 = -subfaultL.depth - np.sin(theta)*subfaultL.depth

        xcenter = 0.5*(xp1 + xp2)
        ycenter = 0.5*(yp1 + yp2)
        fault_width = np.sqrt((xp2-xp1)**2 + (yp2-yp1)**2)

        xcl = xcenter - 0.5*fault_width
        xcr = xcenter + 0.5*fault_width

        self.fault_width = fault_width
        self.fault_depth = -ycenter
        self.xcenter = xcenter
        self.ycenter = ycenter
        self.theta = theta
        self.xcl = xcl
        self.xcr = xcr
        self.xp1 = xp1
        self.xp2 = xp2
        self.yp1 = yp1
        self.yp2 = yp2

    def mapc2p(self,xc,yc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        # constucted signed distance function in computational domain
        ls = numpy.abs(yc - self.ycenter)
        ls = numpy.where(xc < self.xcl, numpy.sqrt((xc-self.xcl)**2 + (yc-self.ycenter)**2), ls)
        ls = numpy.where(xc > self.xcr, numpy.sqrt((xc-self.xcr)**2 + (yc-self.ycenter)**2), ls)

        # define grid that is rotated to line up with fault
        xrot = self.xcenter + numpy.cos(self.theta)*(xc-self.xcenter) + numpy.sin(self.theta)*(yc-self.ycenter)
        yrot = self.ycenter - numpy.sin(self.theta)*(xc-self.xcenter) + numpy.cos(self.theta)*(yc-self.ycenter)

        # Interpolate between rotated grid and cartesian grid near the fault,
        # using cartesian grid far away from fault.
        tol = self.fault_depth
        xp = xc
        yp = yc
        xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xc, xp)
        yp = numpy.where(ls < tol, (tol-ls)/tol*yrot + ls/tol*yc, yp)

        return xp,yp
