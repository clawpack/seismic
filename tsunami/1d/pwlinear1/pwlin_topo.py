
from pylab import *

plot_profile = False

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

if __name__=="__main__":

    from make_nonuniform_grid import make_grid

    dx_min = 3.
    h_min = 5.
    x,B = make_grid(x0,x1,dx_min,h_min,shelf_pwlin,'grid.data')
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


    
    
