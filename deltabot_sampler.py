from __future__ import division
import numpy as np

#deltabot sampler


DELTA_DIAGONAL_ROD=254.2
#// Horizontal offset from middle of printer to smooth rod center.
DELTA_SMOOTH_ROD_OFFSET=175.0
#// Horizontal offset of the universal joints on the end effector.
DELTA_EFFECTOR_OFFSET= 33.0
#// Horizontal offset of the universal joints on the carriages.
DELTA_CARRIAGE_OFFSET=18.0
#// Effective horizontal distance bridged by diagonal push rods.
DELTA_RADIUS= (DELTA_SMOOTH_ROD_OFFSET-DELTA_EFFECTOR_OFFSET-DELTA_CARRIAGE_OFFSET)


MAX_X = 90

R = DELTA_RADIUS
l = DELTA_DIAGONAL_ROD
S = MAX_X
# Tower 3 is at (0, DELTA_RADIUS)
# Max reach between towers, within circle:
rmax1 = min(S, R, l-R)
# Max reach towaards tower
rmax2 = min(S, R, np.sqrt(R**2 + l**2 - R*l))
# but limited by tower! 

nR = 3
dwelltime = 5000 #ms

#########


pointlist = []
def addpoint(x,y):
    pointlist.append((x,y))
    print '%.1f, %.1f'%(x,y)

def addpoints(t0, rvals, nt = 3):
    tvals = t0+np.linspace(0,2*np.pi,nt+1)[0:-1]
    for r in list(rvals):
        for t in list(tvals):
            addpoint(r*np.cos(t),r*np.sin(t))


#complete nR full circles

rvals = np.linspace(0,1.00,nR)[1:]
addpoints(-np.pi/2, rmax1*rvals) #between towers, negative y
addpoints( np.pi/2, rmax2*rvals) #towards towers, positive y
#addpoints(np.pi/2, 100*rvals, 6) 
addpoint(0,0) #origin


import csv

with open('sample.gcode','w') as f:
    f.write(';DELTA_DIAGONAL_ROD= %f\n'%DELTA_DIAGONAL_ROD)
    f.write(';DELTA_RADIUS= %f\n'%DELTA_RADIUS)
    for (x,y) in pointlist:
        f.write('G1 X%06.1f Y%06.1f\nG4 P%.0f\n'%(x,y,dwelltime))

with open('sample.csv','w') as f:
    f.write('DELTA_DIAGONAL_ROD, %f\n'%DELTA_DIAGONAL_ROD)
    f.write('DELTA_RADIUS, %f\n'%DELTA_RADIUS)
    f.writelines(('#%s,\n'%s for s in [
'Data below. Additional comment lines with leading hash-marks may be added',
'At least one comment below parameters. None avove.',
'',
'Only first three columns of data are read. feel free to use other columns',
'Columns are: X, Y, ZHead ( == -Micrometer reading) ',
        ]))
    for (x,y) in pointlist:
            f.write('%6.1f, %6.1f,   0.0\n'%(x,y))
