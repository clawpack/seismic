import clawpack.seismic.dtopotools_horiz_okada as dtopotools
from numpy import arange,sin,pi
reload(dtopotools)

fault = dtopotools.Fault(coordinate_specification='top center')
fault.subfaults = []

width = 50735.0
theta = 0.17
fault_centroid = [25000.0,-19300.0]
slip = 1.0
mu = 3e10
rupture_time = 0.0
rise_time = 1.0
nsubfaults = 50

longitude0 = (fault_centroid[0]-0.5*width)/111.e3
dlongitude = (width/111.e3) / nsubfaults
subfault_width = width/nsubfaults
x = arange(fault_centroid[0]-0.5*width,fault_centroid[0]+0.5*width,subfault_width)
column_list = ['mu','dip','width','depth','slip','rake','strike','length',
                'longitude','latitude','rupture_time','rise_time']

for i in range(nsubfaults):
    subfault = dtopotools.SubFault()
    subfault.mu = mu
    subfault.dip = theta/pi*180.0
    subfault.width = subfault_width
    subfault.depth = -fault_centroid[1] + (x[i]-fault_centroid[0])*sin(theta)
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

fault.write('fault.data',column_list=column_list)
