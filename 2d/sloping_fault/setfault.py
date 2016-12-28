import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from numpy import arange,cos,sin,pi
reload(dtopotools)
from clawpack.geoclaw.data import LAT2METER

fault = dtopotools.Fault(coordinate_specification='top center')
fault.subfaults = []

width = 50000.0
theta = 0.20
fault_top_center = [0.0,-20000.0]
slip = 1.0
mu = 3e10
rupture_time = 0.0
rise_time = 1.0
nsubfaults = 1

longitude0 = fault_top_center[0]/LAT2METER
dlongitude = width*cos(theta)/LAT2METER / nsubfaults
ddepth = width*sin(theta) / nsubfaults
subfault_width = width/nsubfaults

for i in range(nsubfaults):
    subfault = dtopotools.SubFault()
    subfault.mu = mu
    subfault.dip = theta/pi*180.0
    subfault.width = subfault_width
    subfault.depth = -fault_top_center[1] + ddepth*i
    subfault.slip = slip
    subfault.rake = 90
    subfault.strike = 0
    subfault.length = 1000e3
    subfault.longitude = longitude0 + dlongitude*i
    subfault.latitude = 0.
    subfault.coordinate_specification = 'top center'
    subfault.rupture_time = rupture_time
    subfault.rise_time = rise_time

    fault.subfaults.append(subfault)

fault.write('fault.data')
