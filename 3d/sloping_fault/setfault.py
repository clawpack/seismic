import clawpack.seismic.dtopotools_horiz_okada_and_1d as dtopotools
from numpy import arange,sin,pi
reload(dtopotools)

fault = dtopotools.Fault(coordinate_specification='top center')
fault.subfaults = []

width = 50000.0
length = 25000.0
theta = 0.0
fault_centroid = [25000.0,0.0,-20000.0]
slip = 1.0
mu = 3e10
rupture_time = 0.0
rise_time = 1.0
nsubfaults_dip = 40
nsubfaults_strike = 20

longitude0 = (fault_centroid[0]-0.5*width)/111.e3
latitude0 = (fault_centroid[1]-0.5*length)/111.e3
dlongitude = (width/111.e3) / nsubfaults_dip
dlatitude = (length/111.e3) / nsubfaults_strike
subfault_width = width/nsubfaults_dip
subfault_length = length/nsubfaults_strike
x = arange(fault_centroid[0]-0.5*width,fault_centroid[0]+0.5*width,subfault_width)

for i in range(nsubfaults_dip):
    for j in range(nsubfaults_strike):
        subfault = dtopotools.SubFault()
        subfault.mu = mu
        subfault.dip = theta/pi*180.0
        subfault.width = subfault_width
        subfault.depth = -fault_centroid[2] + (x[i]-fault_centroid[0])*sin(theta)
        subfault.slip = slip
        subfault.rake = 90
        subfault.strike = 0
        subfault.length = subfault_length
        subfault.longitude = longitude0 + dlongitude*i
        subfault.latitude = latitude0 + dlatitude*(j+0.5)
        subfault.coordinate_specification = 'top center'
        subfault.rupture_time = rupture_time
        subfault.rise_time = rise_time

        fault.subfaults.append(subfault)

fault.write('fault.data')
