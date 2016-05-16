#!/usr/bin/env python

import numpy as np
from numpy import array as ar

from const_physical import *


### Parameters

boxsize = ar([148.0]*3)
boxcenter = ar([0.0]*3)

net_charge = 0   # AZO:195
#molar_Mg =   5.0e-3 # mol/L (Molarity)
molar_Mg = 0.0
molar_K  = 0.0 # mol/L (Molarity)
#molar_Na = 1.0 # 50.0e-3 # mol/L (Molarity)
#molar_Cl = 1.0 # 50.0e-3 # mol/L (Molarity)
molar_Na = 50.0e-3 # mol/L (Molarity)
molar_Cl = 50.0e-3 # mol/L (Molarity)
filename_base = 'BOX150_NaCl_0050mM'

#############################################

def molar2nion(molar, boxsize):  # boxsize (A)
    #volume = boxsize[0] * 1.0e-10 * boxsize[1] * 1.0e-10 * boxsize[2] * 1.0e-10   (m^3)
    #number = molar (mol/L) * N_AVO (#/mol) * 1.0e3 (L/m^3) * volume (m^3)
    return molar * N_AVO23 * 1.0e-4 * (boxsize[0] * boxsize[1] * boxsize[2])
 
box_lower =  boxcenter + 0.5 * boxsize
box_upper =  boxcenter - 0.5 * boxsize

n_Mg = int(round(molar2nion(molar_Mg,boxsize),0))
n_K  = int(round(molar2nion(molar_K, boxsize),0))
n_Na = int(round(molar2nion(molar_Na,boxsize),0))
n_Cl = (2 * n_Mg + n_K + n_Na - net_charge)

f_log = open('%s.out' % filename_base, 'w')
f_log.write('#Box (x,y,z) = (%f,%f,%f) - (%f,%f,%f)\n' % (box_lower[0],box_lower[1],box_lower[2], box_upper[0],box_upper[1],box_upper[2]))
f_log.write('#Mg: %f mol/L = %f => %i\n' % (molar_Mg, molar2nion(molar_Mg, boxsize), n_Mg))
f_log.write('#K : %f mol/L = %f => %i\n' % (molar_K,  molar2nion(molar_K,  boxsize), n_K))
f_log.write('#Na: %f mol/L = %f => %i\n' % (molar_Na, molar2nion(molar_Na,  boxsize), n_Na))
f_log.write('#Cl: %f mol/L = %f => %i\n' % (molar_Cl, molar2nion(molar_Cl,  boxsize), n_Cl))

xyz_Mg = []
xyz_K  = []
xyz_Na = []
xyz_Cl = []

np.random.seed(111)

def gen_xyz():
    xyz = np.random.random_sample(3)
    sign = np.random.random_sample(3)
    for i,s in enumerate(sign):
        if s < 0.5:
            xyz[i] = -xyz[i]
    xyz = boxcenter + xyz * 0.5 * boxsize 
    return xyz

import energy

ene_pre = 0.0
gap = 2.0

while len(xyz_Mg) < n_Mg:

    tp_Mg = [gen_xyz(),]
    tp_Na = []
    tp_K = []
    tp_Cl = [gen_xyz(), gen_xyz()]

    ene = energy.add_test(xyz_Mg, xyz_K, xyz_Na, xyz_Cl,
                           tp_Mg, tp_K, tp_Na, tp_Cl)

    if ene > gap:
        f_log.write('#      rejection: delta=%f\n' % ene)
    else:
        xyz_Mg.extend(tp_Mg)
        xyz_Cl.extend(tp_Cl)
        ene_pre += ene

    f_log.write('rest: %i   E= %f\n' % (len(xyz_Cl) - n_Cl, ene_pre))


while len(xyz_K) < n_K:

    tp_Mg = []
    tp_K = [gen_xyz(),]
    tp_Na = []
    tp_Cl = [gen_xyz(),]

    ene = energy.add_test(xyz_Mg, xyz_K, xyz_Na, xyz_Cl,
                           tp_Mg, tp_K, tp_Na, tp_Cl)

    if ene > gap:
        f_log.write('#      rejection: delta=%f\n' % ene)
    else:
        xyz_K.extend(tp_K)
        xyz_Cl.extend(tp_Cl)
        ene_pre += ene

    f_log.write('rest: %i   E= %f\n' % (len(xyz_Cl) - n_Cl, ene_pre))


while len(xyz_Na) < n_Na:

    tp_Mg = []
    tp_K = []
    tp_Na = [gen_xyz(),]
    tp_Cl = [gen_xyz(),]

    ene = energy.add_test(xyz_Mg, xyz_K, xyz_Na, xyz_Cl,
                           tp_Mg, tp_K, tp_Na, tp_Cl)

    if ene > gap:
        f_log.write('#      rejection: delta=%f\n' % ene)
    else:
        xyz_Na.extend(tp_Na)
        xyz_Cl.extend(tp_Cl)
        ene_pre += ene

    f_log.write('rest: %i   E= %f\n' % (len(xyz_Cl) - n_Cl, ene_pre))


iatom = 0
f_pdb = open('%s.pdb' % filename_base, 'w')
for xyz in xyz_Mg:
    iatom += 1
    f_pdb.write('ATOM   %4i MG   Mg  B%4i    %8.3f%8.3f%8.3f  1.00  1.00\n' % (iatom,iatom,xyz[0],xyz[1],xyz[2]))
for xyz in xyz_K:
    iatom += 1
    f_pdb.write('ATOM   %4i K    K   B%4i    %8.3f%8.3f%8.3f  1.00  1.00\n' % (iatom,iatom,xyz[0],xyz[1],xyz[2]))
for xyz in xyz_Na:
    iatom += 1
    f_pdb.write('ATOM   %4i Na   Na  B%4i    %8.3f%8.3f%8.3f  1.00  1.00\n' % (iatom,iatom,xyz[0],xyz[1],xyz[2]))
for xyz in xyz_Cl:
    iatom += 1
    f_pdb.write('ATOM   %4i Cl   Cl  B%4i    %8.3f%8.3f%8.3f  1.00  1.00\n' % (iatom,iatom,xyz[0],xyz[1],xyz[2]))
f_pdb.close()

f_log.close()
