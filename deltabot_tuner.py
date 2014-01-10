"""
deltabot_tuner.py: A script for deducing actual deltabot dimensions from micrometer readings.


janus@insignificancegalore.net 
"""



from __future__ import division
from math import sqrt
from collections import namedtuple
import numpy as np
from numpy.linalg import det, norm

## Deltabot basic algebra

Deltabot = namedtuple('Deltabot','tower_origins rod_lenghts')
# tower_origins: tower base coordinates as columns in 3x3 array

def compute_tower_plane_intersection(d1,d2,d3):
    """
    Compute tower-plane intersection from tower separations and the constraints that
    1) circumscribed circle has center in origin 
    2) tower 3 is on positive y axis
    Return as 2x3 array of x and y coordinates for towers
    Tower separations: d_i == |t_{i+1} - t_{i+2}|
    Tower positions: as in Marlin deltabot: 
    1: -x -y
    2: +x -y
    3:  0 +y
    """
    #Radius of circumscribed circle
    R = (d1*d2*d3)/np.sqrt((+d1+d2+d3)*(-d1+d2+d3)*(+d1-d2+d3)*(+d1+d2-d3))
    return np.array((
        (-d2*np.sqrt(1-d2**2/(4*R**2)),  d1*sqrt(1-d1**2/(4*R**2)), 0),
        ((2*R**2 -d2**2)/(2*R), (2*R**2 -d1**2)/(2*R),  R)
        ))

def makebot(d1,d2,d3,z1,z2,z3,l1,l2,l3):
    """
    Construct a Deltabot tuple from tower distances, tower base heights, and arm lengths.
    """
    return Deltabot(
        tower_origins = np.append(compute_tower_plane_intersection(d1,d2,d3), [[z1,z2,z3]], 0),
        rod_lenghts = [float(l) for l in [l1,l2,l3]]
        )

def makebot_ref(d, l):
    """
    Construct a Deltabot tuple for an ideal deltabot with equal tower distances and equal arm lengths.
    """
    return makebot(d,d,d,.0,.0,.0, l,l,l)

def deltabot_head2sled(deltabot, head_position):
    """
    Compute sled position for a certain deltabot geoemetry and head position.

    head_position must be a squeezed vector
    """
    td = deltabot.tower_origins - np.repeat([head_position],3,axis=0).T
    return np.sqrt(
        np.square(deltabot.rod_lenghts) - np.square(td[0:2,:]).sum(0)
            ) - td[2,:]

def norm2(v):
    """For vector: |v|**2"""
    return np.square(v).sum()

def deltabot_sled2head(deltabot, sled_height):
    """
    Compute deltabot head position 
    sled_height:     3-elem vector
    """

    # Compute r: sled position vectors
    sled_position = deltabot.tower_origins.copy()
    sled_position[2,:] = sled_position[2,:]+sled_height

    [r1,r2,r3] = [np.squeeze(r) for r in np.split(sled_position,3,1)]
    [l1,l2,l3] = deltabot.rod_lenghts

    # Compute head position
    D1 = r1-r3
    D2 = r2-r3
    # B ==[D1 D2]
    B = sled_position.dot(np.array(
        ((1,0,-1),(0,1,-1))
        ).T)
    R = B.T.dot(B)
    negt = np.array([[0, 1],[-1, 0]])
    detR = det(R)
    d = np.array([
        norm2(D1) - l1**2 + l3**2,
        norm2(D2) - l2**2 + l3**2
    ])


    s = B.dot(negt.T).dot(R).dot(negt).dot(d) / (2*detR)
    s = r3+s

    hsq4det = 4*detR*np.square(l3) - norm2(B.dot(negt).dot(d))
    normal = np.cross(D1, D2)

    # s is in plane of r1..r3, perpendicularly below head
    # distance to each sled is so that |r_i - s|**2 + h**2 = l_i**2

    return s-(np.sqrt(hsq4det)/(2*detR))*normal 

## Code relating to fit

def roundtrip_error(db0, hp):
    s = deltabot_head2sled(db0, hp)
    return np.linalg.norm(deltabot_sled2head(db0, s) - hp)

from scipy.optimize import minimize
import string

RefPoint = namedtuple('RefPoint', 'xy sled z_meas')
class FitCase:

    def __init__(self, db_ref,  datapoints, db_factory = makebot):
        self.db_ref = db_ref
        self.ref_points = [RefPoint(p[0:2], deltabot_head2sled(db_ref, np.hstack([p[0:2],0])), p[2]) for p in (
            datapoints[i,:] for i in range(datapoints.shape[0]))]
        self.db_factory = db_factory

    def compute_deviations(self, db):
        deviations = [r.z_meas - deltabot_sled2head(db,r.sled)[2] for r in self.ref_points]
        return deviations

    def __call__(self, *args, **kwargs):
        return np.linalg.norm(self.compute_deviations(self.db_factory(*args, **kwargs)))**2

    def do_minimize(self, x0):
        res = minimize(lambda x: self(*x), np.array(x0))
        #print res
        db = self.db_factory(*res.x)
        print 'R: %s'%str(np.sqrt(np.sum(db.tower_origins[0:2,:]**2,0)))
        print 'z: %s'%str(db.tower_origins[2,:])
        print 'l: %s'%str(db.rod_lenghts)
        for p,d in zip(self.ref_points, self.compute_deviations(db)):
            print 'xy: %6.1f %6.1f err: %5.2f  pred: %5.2f'%(p.xy[0],p.xy[1], p.z_meas, d)

def read_FitCase(fit_data):
    """
    Read a data file with header: 
    DELTA_DIAGONAL_ROD, 254.2
    DELTA_RADIUS, xxxx
    #
    """
    parameters = {}
    with open(fit_data) as f:
        #read header
        for line in f:
            if line.startswith('#'):
                break
            k,v = string.split(line,',',maxsplit = 2)[0:2]
            parameters[k.strip()] = float(v)
        db0 = makebot_ref(
            d = parameters['DELTA_RADIUS'] * np.sqrt(3), 
            l = parameters['DELTA_DIAGONAL_ROD'])
        #read data
        datapoints = np.loadtxt(f, delimiter = ',', comments = '#')[:,0:3]
    return FitCase(db0, datapoints)

def main(filename = 'sample.csv'):
    import copy
    fit = read_FitCase(filename)
    d_ref = np.sqrt(3)*norm(fit.db_ref.tower_origins[0:2,0])
    l_ref = fit.db_ref.rod_lenghts[0]

    print 'Z only'
    fit.db_factory = lambda z1,z2,z3: makebot(d_ref,d_ref,d_ref,z1,z2,z3,l_ref,l_ref,l_ref)
    fit.do_minimize(np.array([0,0,0]))

    print 'Common d,l'
    fit.db_factory = lambda d,l,z1,z2,z3: makebot(d,d,d,z1,z2,z3,l,l,l)
    fit.do_minimize(np.array([d_ref, l_ref, 0,0,0]))

    print 'Common l'
    fit2 = copy.copy(fit)
    fit2.db_factory = lambda l,z1,z2,z3: makebot(d_ref, d_ref, d_ref,z1,z2,z3,l,l,l)
    fit2.do_minimize(np.array([l_ref,0,0,0]))

    print 'Individual l'
    fit3 = copy.copy(fit)
    fit3.db_factory = lambda l1,l2,l3,z1,z2,z3: makebot(d_ref, d_ref, d_ref, z1,z2,z3,l1,l2,l3)
    fit3.do_minimize(np.array([l_ref, l_ref, l_ref, 0,0,0]))

    print 'Individual d, l'
    fit4 = copy.copy(fit)
    fit4.db_factory = lambda dd2, dd3, l1,l2,l3,z1,z2,z3: makebot(d_ref, d_ref+dd2, d_ref+dd3, z1,z2,z3,l1,l2,l3)
    fit4.do_minimize(np.array([0,0,l_ref, l_ref, l_ref, 0,0,0]))

if __name__ == '__main__':
    import sys
    if len(sys.argv)==2:
        main(sys.argv[1])
