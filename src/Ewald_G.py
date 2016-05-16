#!/usr/bin/env python

import math
import numpy as np


#L = 150.0
#c = 2.0 * math.pi / L
#print c


'''
h = np.array((0,0,0))
for im in range(0,M):
    for ix in range(0,3):
        print h[0], h[1], h[2], math.sqrt(np.dot(h,h))
        h[ix] += 1
'''
def combsort(x):
    gap = len(x)
    idx = range(gap)
    swapped = True
    while gap > 1 or swapped:
        gap = max(1, int(gap/1.25))
        swapped = False
        for i in range(len(x) - gap):
            j = i+gap
            if x[i] > x[j]:
                x[i], x[j] = x[j], x[i]
                idx[i],idx[j] = idx[j],idx[i]
                swapped = True
    return idx
    

def generate_ReciprocalLattice(M):
    hs = []
    norms = []
    
    for ix in range(-M,M+1):
        for iy in range(-M,M+1):
            for iz in range(-M,M+1):
                if ix==0 and iy==0 and iz==0:
                    continue
                h = np.array((ix,iy,iz))
                norm = math.sqrt(np.dot(h,h))
                hs.append(h)
                norms.append(norm)
    
    idx = combsort(norms)
    
    hs_sorted = []
    for i in range(len(norms)):
        hs_sorted.append(hs[idx[i]])

    return hs_sorted, norms

if __name__ == '__main__':

    hs, norms = generate_ReciprocalLattice(M=20)

    for i in range(len(norms)):
        print '%10i %3i %3i %3i  %12f' % (i, hs[i][0], hs[i][1], hs[i][2], norms[i])
