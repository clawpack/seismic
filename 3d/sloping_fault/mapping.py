import numpy
from pylab import *
import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
reload(dtopotools)
from clawpack.geoclaw.data import LAT2METER

def test(mfault):

    from clawpack.clawutil.data import ClawData

    probdata = ClawData()
    probdata.read('setprob.data',force=True)

    fault = dtopotools.Fault(coordinate_specification='top_center')
    fault.read('fault.data')

    mapping = Mapping(fault)

    domain_depth = probdata.domain_depth
    domain_width = probdata.domain_width
    # num of cells here determined in a similar fashion to that in setrun.py
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
        yp1 = subfaultF.latitude - 0.5*subfaultF.length
        zp1 = -subfaultF.depth
 
        xp2 = subfaultL.longitude*LAT2METER + np.cos(theta)*subfaultL.width
        yp2 = subfaultL.latitude + 0.5*subfaultL.length
        zp2 = -subfaultL.depth - np.sin(theta)*subfaultL.width

        xcenter = 0.5*(xp1 + xp2)
        ycenter = 0.5*(yp1 + yp2)
        zcenter = 0.5*(zp1 + zp2)
        fault_width = np.sqrt((xp2-xp1)**2 + (zp2-zp1)**2)
        fault_length = yp2 - yp1

        xcl = xcenter - 0.5*fault_width
        xcr = xcenter + 0.5*fault_width

        self.fault_width = fault_width
        self.fault_length = fault_length
        self.fault_depth = -zcenter
        self.xcenter = xcenter
        self.ycenter = ycenter
        self.zcenter = zcenter
        self.theta = theta
        self.xcl = xcl
        self.xcr = xcr
        self.xp1 = xp1
        self.xp2 = xp2
        self.zp1 = zp1
        self.zp2 = zp2

	self.slice_xval = None

    def set_slice_xval(self,current_xval):
        self.slice_xval = current_xval

    def mapc2p_xy(self,xc,yc):

        return xc,yc

    def mapc2p_xz(self,xc,zc):
        """
        map computational grid to physical grid that rotates near the fault
        so cell edges match the fault line.  Linear interpolation is used to
        adjust the rotation angle based on distance from fault in computational space.
        The variable tol ensures the physical grid also lines up with a horizontal sea floor
        """

        # constucted signed distance function in computational domain
        ls = numpy.abs(zc - self.zcenter)
        ls = numpy.where(xc < self.xcl, numpy.sqrt((xc-self.xcl)**2 + (zc-self.zcenter)**2), ls)
        ls = numpy.where(xc > self.xcr, numpy.sqrt((xc-self.xcr)**2 + (zc-self.zcenter)**2), ls)

        # define grid that is rotated to line up with fault
        xrot = self.xcenter + numpy.cos(self.theta)*(xc-self.xcenter) + numpy.sin(self.theta)*(zc-self.zcenter)
        zrot = self.zcenter - numpy.sin(self.theta)*(xc-self.xcenter) + numpy.cos(self.theta)*(zc-self.zcenter)

        # Interpolate between rotated grid and cartesian grid near the fault,
        # using cartesian grid far away from fault.
        tol = self.fault_depth
        xp = xc
        zp = zc
        xp = numpy.where(ls < tol, (tol-ls)/tol*xrot + ls/tol*xc, xp)
        zp = numpy.where(ls < tol, (tol-ls)/tol*zrot + ls/tol*zc, zp)

        return xp,zp

    def mapc2p_yz(self,yc,zc):

        xc = self.slice_xval

        # constucted signed distance function in computational domain
        if (xc < self.xcl):
            ls = numpy.sqrt((xc-self.xcl)**2 + (zc-self.zcenter)**2)
        elif (xc > self.xcr):
            ls = numpy.sqrt((xc-self.xcr)**2 + (zc-self.zcenter)**2)
        else:
            ls = numpy.abs(zc - self.zcenter)

        zrot = self.zcenter - numpy.sin(self.theta)*(xc-self.xcenter) + numpy.cos(self.theta)*(zc-self.zcenter)

        # Interpolate between rotated grid and cartesian grid near the fault,
        # using cartesian grid far away from fault.
        tol = self.fault_depth
        yp = yc
        zp = numpy.where(ls < tol, (tol-ls)/tol*zrot + ls/tol*zc, zc)

        return yp,zp
