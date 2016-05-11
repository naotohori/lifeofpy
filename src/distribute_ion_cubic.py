#!/usr/bin/env python

import numpy as np
from numpy import array as ar

from const_physical import *

#boxsize = ar([350.0]*3)
boxsize = ar([150.0]*3)
boxcenter = ar([0.0]*3)

box_lower =  boxcenter + 0.5 * boxsize
box_upper =  boxcenter - 0.5 * boxsize

print '#Box', box_lower, box_upper

net_charge = 0   # AZO:195
#molar_Mg =   5.0e-3 # mol/L (Molarity)
molar_Mg = 0.0
molar_K  = 0.0 # mol/L (Molarity)
molar_Na = 1.0 # 50.0e-3 # mol/L (Molarity)
molar_Cl = 1.0 # 50.0e-3 # mol/L (Molarity)

def molar2nion(molar, boxsize):  # boxsize (A)
    #volume = boxsize[0] * 1.0e-10 * boxsize[1] * 1.0e-10 * boxsize[2] * 1.0e-10   (m^3)
    #number = molar (mol/L) * N_AVO (#/mol) * 1.0e3 (L/m^3) * volume (m^3)
    return molar * N_AVO23 * 1.0e-4 * (boxsize[0] * boxsize[1] * boxsize[2])


#n_Mg = int(round(molar2nion(molar_Mg,boxsize),0))
n_Mg = 0
#n_K  = int(round(molar2nion(molar_K, boxsize),0))
n_K  = 0
n_Na = int(round(molar2nion(molar_Na,boxsize),0))
n_Cl = (2 * n_Mg + n_K + n_Na - net_charge)
print '#Mg', molar_Mg, molar2nion(molar_Mg, boxsize), n_Mg
print '#K ', molar_K,  molar2nion(molar_K,  boxsize), n_K
print '#Na', molar_Na, molar2nion(molar_Na, boxsize), n_Na
print '#Cl', molar_Cl, molar2nion(molar_Cl, boxsize), n_Cl

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

ene_pre = 1.0
gap = 2.0

while len(xyz_Mg) < n_Mg:

    xyz_Mg.append(gen_xyz())
    xyz_Cl.append(gen_xyz())
    xyz_Cl.append(gen_xyz())

    ene = energy.sumup(xyz_Mg, xyz_K, xyz_Na, xyz_Cl)
    print('rest: %i   E= %f' % (len(xyz_Cl) - n_Cl, ene))

    if ene > ene_pre+gap:
        xyz_Mg.pop()
        xyz_Cl.pop()
        xyz_Cl.pop()
        print('#      rejected')
    else:
        ene_pre = ene


while len(xyz_K) < n_K:

    xyz_K.append(gen_xyz())
    xyz_Cl.append(gen_xyz())

    ene = energy.sumup(xyz_Mg, xyz_K, xyz_Na, xyz_Cl)
    print('rest: %i   E= %f' % (len(xyz_Cl) - n_Cl, ene))

    if ene > ene_pre+gap:
        xyz_K.pop()
        xyz_Cl.pop()
        print('#      rejected')
    else:
        ene_pre = ene


while len(xyz_Na) < n_Na:

    xyz_Na.append(gen_xyz())
    xyz_Cl.append(gen_xyz())

    ene = energy.sumup(xyz_Mg, xyz_K, xyz_Na, xyz_Cl)
    print('rest: %i   E= %f' % (len(xyz_Cl) - n_Cl, ene))

    if ene > ene_pre+gap:
        xyz_Na.pop()
        xyz_Cl.pop()
        print('#      rejected')
    else:
        ene_pre = ene

print '#Mg', len(xyz_Mg)
print '#K ', len(xyz_K)
print '#Na', len(xyz_Na)
print '#Cl', len(xyz_Cl)

iatom = 0
#f_pdb = open('ion.pdb', 'w')
#f_pdb = open('BOX150_NaCl_0050mM.pdb', 'w')
f_pdb = open('BOX150_NaCl_1000mM.pdb', 'w')
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
