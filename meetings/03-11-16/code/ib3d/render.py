#!/usr/bin/env python
import sys
import os

path = "./op_sphere/"
pv = 'python vtk2pov.py'
if not os.path.exists("./rendered/"):
    os.makedirs("./rendered/")

if __name__ == '__main__':
    num_frame = len([f for f in os.listdir(path)
                if os.path.isfile(os.path.join(path, f))])
    for frame in xrange(num_frame):
        Xmemb = path + "%03d.vtk"%frame
        os.system('{0} {1}'.format(pv, Xmemb))
        os.rename("%03d.inc"%frame, 'Xmemb.inc')
        os.system('povray seafloor.pov -D')
        # High quality: os.system('povray +A0.01 +R9 +Q9 -J seafloor.pov -D')
        os.rename('seafloor.png', './rendered/%03d.png'%frame)
