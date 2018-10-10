# -*- coding: utf-8 -*-
"""Test Program for spktype01.

This program computes differences of the result of spktype01 from the 
ephemeris provided by JPL HORIZONS system with three celestial bodies: 
dwarf planet Ceres, and asteroid Vesta.

This program reads two groups of data files.
First group is SPK data files; Ceres_Vesta.bsp.  this file was 
created by HORIZONS system upon user requests through telnet.
Second group is ephemeris data files; Ceres_jpl.csv, and Vesta_jpl.csv.  
Each file contains sun centered positions and velocities of respective 
celestial body.  These data are copied from outputs of HORIZONS system.

This program computes 'ranges' for position and velocity, and compare them to
limits, dposlim and dvellim.  The limits are arbitrary, and physical meanings
are unclear.

@author: Shushi Uetsuki (whiskie14142)
"""

from spktype01 import SPKType01
import csv
import numpy as np

bspfiles = 'Ceres_Vesta.bsp'
csvfiles = ['Ceres_jpl.csv', 'Vesta_jpl.csv']
testnames = ['Ceres', 'Vesta']
bodyids = [2000001, 2000004]

dposlim = 1.0       # position difference limit = 1.0 kilometer
dvellim = 1e-7      # velocity difference limit = 0.0001 m/s

kernel = SPKType01.open(bspfiles)
print(kernel)

for testcase in range(2):
    csvfile = open(csvfiles[testcase])
    print()
    print(testnames[testcase])
    print()
    
    err = 0
    count = 0
    for row in csv.reader(csvfile):
        count += 1
        refjd = float(row[0])
        refpos = np.array([float(row[2]), float(row[3]), float(row[4])])
        refvel = np.array([float(row[5]), float(row[6]), float(row[7])])
        
        spkpos, spkvel = kernel.compute_type01(10, bodyids[testcase], refjd)
        
        dpos = refpos - spkpos
        dvel = refvel - spkvel
        prange = np.sqrt(np.dot(dpos, dpos))
        vrange = np.sqrt(np.dot(dvel, dvel))
        if prange > dposlim or vrange > dvellim:
            print(count, refjd, prange, vrange, ' *')
            err += 1
        else:
            print(count, refjd, prange, vrange)
    print('Checked count = ', count, ',  Error count = ', err)
    csvfile.close()
kernel.close()
