
from pylab import *
from scipy.interpolate import interp1d


plot_profile = True

grav = 9.81
mx = 10000   # number of grid cells

x0 = -150e3           # left boundary (meters)
x0_slope = -65e3      # start of slope
x0_shelf = -45e3      # start of shelf
x0_beach = -5e3       # start of beach
x0_shore = 0.         # initial shoreline
x1 = x0_shore + 2e3   # right boundary

z0_ocean = -3000.     # depth of ocean
z0_shelf = -100.      # depth at x0_shelf
z0_beach = -100.       # depth at x0_beach
z0_shore = 0.         # depth at x0_shore

## Used by sloping_fault code to define seafloor so topo matches
def get_seafloor_parameters():
    return x0, x0_slope, x0_shelf, x0_beach, x0_shore, x1

if x0_beach != x0_shelf:
    slope_of_shelf = (z0_beach - z0_shelf) / (x0_beach - x0_shelf)
else:
    slope_of_shelf = 1e9

if x0_slope != x0_shelf:
    slope_of_slope = (z0_ocean - z0_shelf) / (x0_slope - x0_shelf)
else:
    slope_of_slope = 1e9

slope_of_beach = (z0_beach - z0_shore) / (x0_beach - x0_shore)
print "Slope of shelf = ",slope_of_shelf
print "Slope of beach = ",slope_of_beach

def shelf_pwlin(r):
    """
    Ocean followed by continental slope, continental shelf, and beach.
    The ocean is flat, the slope, shelf, and beach are linear.
    """
    z = z0_shore + slope_of_beach * r   # beach
    z = where(r<x0_beach, z0_shelf + slope_of_shelf*(r-x0_shelf), z)
    z = where(r<x0_shelf, z0_ocean + slope_of_slope*(r-x0_slope), z)
    z = where(r<x0_slope, z0_ocean, z)
    return z

if plot_profile:
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

# Make the nonuniform grid:

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
    z = shelf_pwlin(x)
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

dxmin = diff(xp).min()
dxmax = diff(xp).max()
print "Maximum cell width is %7.2f m, minimum cell width is %7.2f m" \
            % (dxmax,dxmin)

f = open('grid.data','w')
f.write('%10i \n' % mx)
for i in range(mx+1):
    f.write('%15.4f %15.4f\n' % (xp[i],z[i]))
f.close()
