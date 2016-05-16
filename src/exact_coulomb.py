#!/usr/bin/env python

import numpy as np
import energy
import Ewald_G
from const_physical import *

TK = 310.0
MM_A = 87.740
MM_B = -0.4008
MM_C = 9.398e-4
MM_D = -1.410e-6
Tc = TK - 273.15
ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

L = 150.0

def calc_coulomb(xyz,charge,hvec,flg_self):
    
    ene = 0.0 
    nmp = len(xyz)

    # Within the same box
    #if hvec[0] == 0 and hvec[1] == 0 and hvec[2] == 0:
    if flg_self:
        for i in range(nmp):
            for j in range(i+1,nmp):
                ene += energy.ele_coulomb(xyz[i],xyz[j], charge[i]*charge[j],-1.0)
    
    else:
        for j in range(nmp):
            xyzj = xyz[j] + hvec * L
            qj = charge[j]
            for i in range(nmp):
                ene += energy.ele_coulomb(xyz[i], xyzj, qj*charge[i], -1.0)
    
    return ene

if __name__ == '__main__':
    xyz = []
    charge = []
    
    for i,l in enumerate(open('BOX150_NaCl_0050mM.pdb','r')):
        xyz.append(np.array( (float(l[30:38]),float(l[38:46]),float(l[46:54])) ))
        if i < 98:
            charge.append(+1.0)
        else:
            charge.append(-1.0)
    
    
    coef = JOUL2KCAL_MOL * 1.0e10 * ELE*ELE / (4.0*F_PI*EPSI_0*ek) 
    print 'coef=',coef

    hvecs, norms = Ewald_G.generate_ReciprocalLattice(50)

    #hmax = 50
    #hmax2 = hmax * hmax
    #hvecs = []
    #for iL in range(-hmax,hmax):
    #    for jL in range(-hmax,hmax):
    #        for kL in range(-hmax,hmax):
    #            if iL==0 and jL==0 and kL==0:
    #                continue
    #            norm = iL**2 + jL**2 + kL**2
    #            if norm <= hmax2:
    #                hvecs.append(np.array((iL,jL,kL)))
    #print len(hvecs)


    ene = coef * calc_coulomb(xyz,charge,(0.0,0.0,0.0),True)
    ene_total = ene
    print ('   #%8.3f %4i %4i %4i %10.4f %10.4f' % (0.0,0,0,0,ene,ene_total))
    print ('%8.3f %10.4f' % (0.0,ene_total))
    norm_pre = 0.0
    
    for hvec,norm in zip(hvecs,norms):
        ene = coef * calc_coulomb(xyz,charge,hvec,False)
        ene_total += ene
        print ('   #%8.3f %4i %4i %4i %10.4f %10.4f' % (norm,hvec[0],hvec[1],hvec[2],ene,ene_total))

        if norm != norm_pre:
            print ('%8.3f %10.4f' % (norm,ene_total))
            norm_pre = norm
