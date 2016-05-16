#!/usr/bin/env python

import numpy as np
from math import sqrt
from param import *

def ele_coulomb(xyz1, xyz2, coef, cutoff):
    v21 = xyz2 - xyz1
    dist2 = np.dot(v21, v21)
    if cutoff > 0 and dist2 > (cutoff ** 2):
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
    d = exv_adjust / dist
    d2 = d  * d
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


def add_test(xyz_Mg, xyz_K, xyz_Na, xyz_Cl,
              tp_Mg, tp_K, tp_Na, tp_Cl):
    ene = 0.0

    nMg = len(xyz_Mg)
    nK  = len(xyz_K )
    nNa = len(xyz_Na)
    nCl = len(xyz_Cl)

    ntpMg = len(tp_Mg)
    ntpK  = len(tp_K )
    ntpNa = len(tp_Na)
    ntpCl = len(tp_Cl)

    for i in range(ntpMg):
        for j in range(nMg):
            ene += ele_coulomb(tp_Mg[i], xyz_Mg[j], coef_ele_Mg_Mg)
            e = exv_dt15(tp_Mg[i], xyz_Mg[j], coef_exv_Mg_Mg, dist_exv_Mg_Mg)
            if e > 999999:
                return e
            ene += e

        for j in range(nK):
            ene += ele_coulomb(tp_Mg[i], xyz_K[j], coef_ele_Mg_K)
            e = exv_dt15(tp_Mg[i], xyz_K[j], coef_exv_Mg_K, dist_exv_Mg_K)
            if e > 999999:
                return e
            ene += e

        for j in range(nNa):
            ene += ele_coulomb(tp_Mg[i], xyz_Na[j], coef_ele_Mg_Na)
            e = exv_dt15(tp_Mg[i], xyz_Na[j], coef_exv_Mg_Na, dist_exv_Mg_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(nCl):
            ene += ele_coulomb(tp_Mg[i], xyz_Cl[j], coef_ele_Mg_Cl)
            e = exv_dt15(tp_Mg[i], xyz_Cl[j], coef_exv_Mg_Cl, dist_exv_Mg_Cl)
            if e > 999999:
                return e
            ene += e

        for j in range(i+1, ntpMg):
            ene += ele_coulomb(tp_Mg[i], tp_Mg[j], coef_ele_Mg_Mg)
            e = exv_dt15(tp_Mg[i], tp_Mg[j], coef_exv_Mg_Mg, dist_exv_Mg_Mg)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpK):
            ene += ele_coulomb(tp_Mg[i], tp_K[j], coef_ele_Mg_K)
            e = exv_dt15(tp_Mg[i], tp_K[j], coef_exv_Mg_K, dist_exv_Mg_K)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpNa):
            ene += ele_coulomb(tp_Mg[i], tp_Na[j], coef_ele_Mg_Na)
            e = exv_dt15(tp_Mg[i], tp_Na[j], coef_exv_Mg_Na, dist_exv_Mg_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpCl):
            ene += ele_coulomb(tp_Mg[i], tp_Cl[j], coef_ele_Mg_Cl)
            e = exv_dt15(tp_Mg[i], tp_Cl[j], coef_exv_Mg_Cl, dist_exv_Mg_Cl)
            if e > 999999:
                return e
            ene += e

    for i in range(ntpK):

        for j in range(nK):
            ene += ele_coulomb(tp_K[i], xyz_K[j], coef_ele_K_K)
            e = exv_dt15(tp_K[i], xyz_K[j], coef_exv_K_K, dist_exv_K_K)
            if e > 999999:
                return e
            ene += e

        for j in range(nNa):
            ene += ele_coulomb(tp_K[i], xyz_Na[j], coef_ele_K_Na)
            e = exv_dt15(tp_K[i], xyz_Na[j], coef_exv_K_Na, dist_exv_K_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(nCl):
            ene += ele_coulomb(tp_K[i], xyz_Cl[j], coef_ele_K_Cl)
            e = exv_dt15(tp_K[i], xyz_Cl[j], coef_exv_K_Cl, dist_exv_K_Cl)
            if e > 999999:
                return e
            ene += e

        for j in range(i+1,ntpK):
            ene += ele_coulomb(tp_K[i], tp_K[j], coef_ele_K_K)
            e = exv_dt15(tp_K[i], tp_K[j], coef_exv_K_K, dist_exv_K_K)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpNa):
            ene += ele_coulomb(tp_K[i], tp_Na[j], coef_ele_K_Na)
            e = exv_dt15(tp_K[i], tp_Na[j], coef_exv_K_Na, dist_exv_K_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpCl):
            ene += ele_coulomb(tp_K[i], tp_Cl[j], coef_ele_K_Cl)
            e = exv_dt15(tp_K[i], tp_Cl[j], coef_exv_K_Cl, dist_exv_K_Cl)
            if e > 999999:
                return e
            ene += e

    for i in range(ntpNa):

        for j in range(nNa):
            ene += ele_coulomb(tp_Na[i], xyz_Na[j], coef_ele_Na_Na)
            e = exv_dt15(tp_Na[i], xyz_Na[j], coef_exv_Na_Na, dist_exv_Na_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(nCl):
            ene += ele_coulomb(tp_Na[i], xyz_Cl[j], coef_ele_Na_Cl)
            e = exv_dt15(tp_Na[i], xyz_Cl[j], coef_exv_Na_Cl, dist_exv_Na_Cl)
            if e > 999999:
                return e
            ene += e

        for j in range(i+1,ntpNa):
            ene += ele_coulomb(tp_Na[i], tp_Na[j], coef_ele_Na_Na)
            e = exv_dt15(tp_Na[i], tp_Na[j], coef_exv_Na_Na, dist_exv_Na_Na)
            if e > 999999:
                return e
            ene += e

        for j in range(ntpCl):
            ene += ele_coulomb(tp_Na[i], tp_Cl[j], coef_ele_Na_Cl)
            e = exv_dt15(tp_Na[i], tp_Cl[j], coef_exv_Na_Cl, dist_exv_Na_Cl)
            if e > 999999:
                return e
            ene += e

    for i in range(ntpCl):
        for j in range(nNa):
            ene += ele_coulomb(tp_Cl[i], xyz_Na[j], coef_ele_Na_Cl)
            e = exv_dt15(tp_Cl[i], xyz_Na[j], coef_exv_Na_Cl, dist_exv_Na_Cl)
            if e > 999999:
                return e
            ene += e

        for j in range(nCl):
            ene += ele_coulomb(tp_Cl[i], xyz_Cl[j], coef_ele_Cl_Cl)
            e = exv_dt15(tp_Cl[i], xyz_Cl[j], coef_exv_Cl_Cl, dist_exv_Cl_Cl)
            if e > 999999:
                return e
            ene += e

        for j in range(i+1,ntpCl):
            ene += ele_coulomb(tp_Cl[i], tp_Cl[j], coef_ele_Cl_Cl)
            e = exv_dt15(tp_Cl[i], tp_Cl[j], coef_exv_Cl_Cl, dist_exv_Cl_Cl)
            if e > 999999:
                return e
            ene += e

    return ene
