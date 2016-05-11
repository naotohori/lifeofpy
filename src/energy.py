#!/usr/bin/env python

import numpy as np
from math import sqrt
from param import *

def ele_coulomb(xyz1, xyz2, coef):
    v21 = xyz2 - xyz1
    dist2 = np.dot(v21, v21)
    if dist2 > (999.0 ** 2):
        return 0.0
    dist = sqrt(dist2)
    return coef / dist


def exv_dt15(xyz1, xyz2, coef, d_rad):
    v21 = xyz2 - xyz1
    dist = sqrt(np.dot(v21,v21))
    if dist > d_rad:
        return 0.0
    dist = dist + exv_adjust - d_rad
    if dist < exv_inf:
        return 99999999
    d2 = exv_adjust / dist
    d4 = d2 * d2
    d6 = d2 * d4
    d12= d6 * d6
    return coef * (d12 - 2*d6 + 1.0)


def sumup(xyz_Mg, xyz_K, xyz_Na, xyz_Cl):
    ene = 0.0

    nMg = len(xyz_Mg)
    nK  = len(xyz_K )
    nNa = len(xyz_Na)
    nCl = len(xyz_Cl)

    for i in range(nMg-1):
        for j in range(i+1,nMg):
            ene += ele_coulomb(xyz_Mg[i], xyz_Mg[j], coef_ele_Mg_Mg)
            e = exv_dt15(xyz_Mg[i], xyz_Mg[j], coef_exv_Mg_Mg, dist_exv_Mg_Mg)
            if e > 999999:
                return e
            ene += e


    for i in range(nK-1):
        for j in range(i+1,nK):
            ene += ele_coulomb(xyz_K[i], xyz_K[j], coef_ele_K_K)
            e = exv_dt15(xyz_K[i], xyz_K[j], coef_exv_K_K, dist_exv_K_K)
            if e > 999999:
                return e
            ene += e

    for i in range(nNa-1):
        for j in range(i+1,nNa):
            ene += ele_coulomb(xyz_Na[i], xyz_Na[j], coef_ele_Na_Na)
            e = exv_dt15(xyz_Na[i], xyz_Na[j], coef_exv_Na_Na, dist_exv_Na_Na)
            if e > 999999:
                return e
            ene += e

    for i in range(nCl-1):
        for j in range(i+1,nCl):
            ene += ele_coulomb(xyz_Cl[i], xyz_Cl[j], coef_ele_Cl_Cl)
            e = exv_dt15(xyz_Cl[i], xyz_Cl[j], coef_exv_Cl_Cl, dist_exv_Cl_Cl)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_Mg:
        for xyz2 in xyz_K:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_Mg_K)
            e = exv_dt15(xyz1, xyz2, coef_exv_Mg_K, dist_exv_Mg_K)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_Mg:
        for xyz2 in xyz_Na:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_Mg_Na)
            e = exv_dt15(xyz1, xyz2, coef_exv_Mg_Na, dist_exv_Mg_Na)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_Mg:
        for xyz2 in xyz_Cl:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_Mg_Cl)
            e = exv_dt15(xyz1, xyz2, coef_exv_Mg_Cl, dist_exv_Mg_Cl)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_K:
        for xyz2 in xyz_Na:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_K_Na)
            e = exv_dt15(xyz1, xyz2, coef_exv_K_Na, dist_exv_K_Na)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_K:
        for xyz2 in xyz_Cl:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_K_Cl)
            e = exv_dt15(xyz1, xyz2, coef_exv_K_Cl, dist_exv_K_Cl)
            if e > 999999:
                return e
            ene += e

    for xyz1 in xyz_Na:
        for xyz2 in xyz_Cl:
            ene += ele_coulomb(xyz1, xyz2, coef_ele_Na_Cl)
            e = exv_dt15(xyz1, xyz2, coef_exv_Na_Cl, dist_exv_Na_Cl)
            if e > 999999:
                return e
            ene += e
    
    return ene
