#!/usr/bin/env python

from math import sqrt 
from const_physical import ELE, JOUL2KCAL_MOL, F_PI, EPSI_0

exv_coef = 1.0
exv_adjust = 1.5852
exv_dist = 3.2
exv_inf  = 0.2

#Mg2  0.7926  0.894700
#Ca2  1.7131  0.459789
#Cl   1.9480  0.265000
#K    2.6580  0.000328
#Na   1.8680  0.002770
exv_Mg_R = 0.7926
exv_K_R  = 2.6580
exv_Na_R = 1.8680
exv_Cl_R = 1.9480

exv_Mg_eps = 0.894700
exv_K_eps  = 0.000328
exv_Na_eps = 0.002770
exv_Cl_eps = 0.265000

charge_Mg = +2.0
charge_K  = +1.0
charge_Na = +1.0
charge_Cl = -1.0

TK = 310.0
MM_A = 87.740
MM_B = -0.4008
MM_C = 9.398e-4
MM_D = -1.410e-6
Tc = TK - 273.15
ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

coef_ele_Mg_Mg = charge_Mg * charge_Mg * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_K_K   = charge_K  * charge_K  * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Na_Na = charge_Na * charge_Na * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Cl_Cl = charge_Cl * charge_Cl * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Mg_K  = charge_Mg * charge_K  * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Mg_Na = charge_Mg * charge_Na * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Mg_Cl = charge_Mg * charge_Cl * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_K_Na  = charge_K  * charge_Na * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_K_Cl  = charge_K  * charge_Cl * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)
coef_ele_Na_Cl = charge_Na * charge_Cl * JOUL2KCAL_MOL * 1.0e10 * ELE * ELE / (4.0 * F_PI * EPSI_0 * ek)

coef_exv_Mg_Mg = sqrt(exv_Mg_eps * exv_Mg_eps)
coef_exv_K_K   = sqrt(exv_K_eps  * exv_K_eps)
coef_exv_Na_Na = sqrt(exv_Na_eps * exv_Na_eps)
coef_exv_Cl_Cl = sqrt(exv_Cl_eps * exv_Cl_eps)
coef_exv_Mg_K  = sqrt(exv_Mg_eps * exv_K_eps)
coef_exv_Mg_Na = sqrt(exv_Mg_eps * exv_Na_eps)
coef_exv_Mg_Cl = sqrt(exv_Mg_eps * exv_Cl_eps)
coef_exv_K_Na  = sqrt(exv_K_eps  * exv_Na_eps)
coef_exv_K_Cl  = sqrt(exv_K_eps  * exv_Cl_eps)
coef_exv_Na_Cl = sqrt(exv_Na_eps * exv_Cl_eps)

dist_exv_Mg_Mg = exv_Mg_R * 2
dist_exv_K_K   = exv_K_R  * 2
dist_exv_Na_Na = exv_Na_R * 2
dist_exv_Cl_Cl = exv_Cl_R * 2
dist_exv_Mg_K  = exv_Mg_R + exv_K_R
dist_exv_Mg_Na = exv_Mg_R + exv_Na_R
dist_exv_Mg_Cl = exv_Mg_R + exv_Cl_R
dist_exv_K_Na  = exv_K_R  + exv_Na_R
dist_exv_K_Cl  = exv_K_R  + exv_Cl_R
dist_exv_Na_Cl = exv_Na_R  + exv_Cl_R
