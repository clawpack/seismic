
"""
Code to create a nonuniform grid for 1d geoclaw modeling with the property
that the cell width is greater in deep water and smaller in shallow water,
so that the Courant number is roughly constant up to some depth off shore
and then the grid is uniform to the shoreline and beyond.  
"""

from pylab import *
from scipy.interpolate import interp1d

grav = 9.81   # wave speed is approx sqrt(grav * h) in water of depth h


if 0:
    r = linspace(x0,x1,1000)
    s = shelf_pwlin(r)
    eta = where(s<0,0,s)

    figure(13,figsize=(12,5))
    clf()
    fill_between(r/1e3,s,eta,color='b')
    plot(r/1e3,s,'g')
    xlim(x0/1e3, x1/1e3)
    ylim(z0_ocean*1.1, 200)
    xlabel('kilometers')
    ylabel('meters')
    title('shelf and beach profile')
    fname = 'profile.png'
    savefig(fname)
    print "Created ",fname

def make_grid(x1, x2, dx_min, h_min, topo_func, fname=None):
    """
    :Input:
        - `x1, x2, mx:` Create a grid with mx cells between x1 and x2
        - `dx_onshore:` Minimum cell width used on dry land and nearshore
        - `topo_func:` Function the returns the topo B(x) for any x
    :Output:
        - `xc`: `mx+1` equally spaced computational points
        - `xp`: `mx+1` cell edges with unequal spacing
    """

    sea_level = 0.
    dt = dx_min / sqrt(grav*h_min)
    x = [x1]
    for j in range(100000):
        xj = x[-1]
        Bj = topo_func(xj)
        hj = max(sea_level - Bj, 0.)
        dxj = max(dx_min, dt*sqrt(grav*hj))
        xjp1 = xj + dxj

        # check sqrt(gh) at right edge of provisional cell in case that's
        # smaller:

        Bjp1 = topo_func(xjp1)
        hjp1 = max(sea_level - Bjp1, 0.)
        dxj = max(dxj, dt*sqrt(grav*hjp1))
        xjp1 = xj + dxj

        if xjp1 <= x2:
            x.append(xj + dxj)
        else:
            # increase size of final cell to reach right boundary:
            x[-1] = x2
            break
            
    x = array(x)
    mx = len(x) - 1 # number of grid cells
    print "Created grid with %i cells, for dt <= %6.3f" % (mx,dt)
    dxmin = diff(x).min()
    dxmax = diff(x).max()
    print "Maximum cell width is %7.2f m, minimum cell width is %7.2f m" \
                % (dxmax,dxmin)
    B = topo_func(x)

    if fname is not None:
        f = open(fname,'w')
        f.write('%10i \n' % mx)
        for i in range(mx+1):
            f.write('%15.4f %15.4f\n' % (x[i],B[i]))
        f.close()
        print "Created %s" % fname
        
    return x,B
    
    
if 0:

    # xc will be uniformly spaced computational grid
    # xp will be physical grid of cell edges from x0 to x1 as defined above

    # The grid will be nonuniform and chosen so that
    #   c(xp[i]) / (xp[i+1] - xp[i]) \approx constant
    # where c(x) = sqrt(grav * h(x)) is the wave speed
    # so that the Courant number is roughly constant everywhere.

    # But near shore the depth goes to zero, so set a minimum depth
    # that is used in calculating c, and then the Courant number will be
    # roughly constant in deeper water but the cells will be uniform in
    # shallower water and the Courant number will decrease as shore is approached.

    hmin = 50.  # minimum depth to use in computing c
    cmin = sqrt(grav*hmin)

    def c(x):
        z = topo_func(x)
        h = where(-z > hmin, -z, hmin)
        c = sqrt(grav*h)
        return c

    xunif = linspace(x0, x1, 2*mx)
    cunif = c(xunif)
    csum = cumsum(1./cunif)
    csum = csum - csum[0]

    csum = csum / csum[-1]
    cinv = interp1d(csum, xunif)

    xc = linspace(0, 1, mx+1)   # computational grid
    xp = cinv(xc)
    z = shelf_pwlin(xp)


    f = open('grid.data','w')
    f.write('%10i \n' % mx)
    for i in range(mx+1):
        f.write('%15.4f %15.4f\n' % (xp[i],z[i]))
    f.close()
