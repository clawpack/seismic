#!/usr/bin/env python
r"""Slice data class for 3d output in slices
"""

from clawpack.clawutil.data import ClawData


class SliceData(ClawData):

    def __init__(self):
        super(SliceData,self).__init__()

        self.add_attribute('slices',[])

    def add(self,point,normal):
        new_slice = ClawSlice(point,normal)
        self.slices.append(new_slice)

    def write(self,out_file='slices.data',data_source='setrun.py'):
        r"""Write out slice information data file."""

        # Write out slice data file
        self.open_data_file(out_file,data_source)
        self.data_write(name='nslices', value=len(self.slices))
        for current_slice in self.slices:
            self.data_write(name='point',value=current_slice.point)
            self.data_write(name='normal',value=current_slice.normal)
        self.close_data_file()

    def read(Self,data_path="./",filename='slices.data'):
        r"""Read slice data file"""
        path = os.path.join(data_path, file_name)
        slice_file = open(path,'r')

        # Read past comments and blank lines
        header_lines = 0
        ignore_lines = True
        while ignore_lines:
            line = slice_file.readline()
            if line[0] == "#" or len(line.strip()) == 0:
                header_lines += 1
            else:
                break

        # Read number of slices
        num_slices = int(line.split()[0])

        # Read slices
        for n in xrange(num_slices):
            line = slice_file.readline().split()
            point = [float(a) for a in line[3:6]]
            normal = [float(a) for a in line[0:3]]
            self.add(point,normal)

        slice_file.close()

class ClawSlice():
    r"""An output plane defined by a point and normal direction"""

    def __init__(self,point,normal):
        self.point = point
        self.normal = normal
