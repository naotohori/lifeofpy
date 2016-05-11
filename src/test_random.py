#!/usr/bin/env python

import numpy as np

np.random.seed(100)

for it in range(1000):
    xyz = np.random.random_sample(3)
    sign = np.random.random_sample(3)
    
    for i,s in enumerate(sign):
        if s < 0.5:
            xyz[i] = -xyz[i]

    #xyz = xyz * 350 / 2. + 0
    xyz = xyz * (185 - (-165)) / 2. + 10
    
    print xyz[0], xyz[1], xyz[2]
    
