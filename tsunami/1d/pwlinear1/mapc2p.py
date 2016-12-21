
import numpy as np
d = np.loadtxt('grid.data', skiprows=1)
print 'Read grid from grid.data, %i grid values' % d.shape[0]

def mapc2p(xc):
    xp_edges = d[:,0]
    xp = 0.5*(xp_edges[:-1] + xp_edges[1:])
    if len(xp) != len(xc):
        print '*** Length of xp from grid.data does not match xc'
        raise ValueError('*** Length of xp from grid.data (%i) does not match xc (%i)' \
                % (len(xp),len(xc)))
    return xp
