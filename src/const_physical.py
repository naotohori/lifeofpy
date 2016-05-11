#module const_physical
#  use const_maxsize

SDIM = 3 # # of space dimension

F_PI    = 3.14159265358979323846264338e0 # Circular constant (Pi)
F_2PI   = F_PI + F_PI

EPSI_0  = 8.854187817e-12   # Vacuum permittivity [F/m]
K_BOLTZ = 1.380650424e-23   # Boltzmann constant [J/K]
BOLTZC  = 1.986231313e-3    # Boltzmann constant [kcal/mol/K]
##  real(PREC), parameter :: BOLTZC  = 1.98e-3    # Boltzmann constant [kcal/mol/K]
N_AVO   = 6.0221417930e23   # Avogadro constant
N_AVO23 = 6.0221417930      # Avogadro constant
ELE     = 1.60217648740e-19 # Elementary charge [C]
JOUL2KCAL_MOL = 1.43862e20  # (J -> kcal/mol)
JOUL2KCAL = 4186.8   # (J -> kcal)

DE_MAX  = 20.0e0 # limit value of force

  # judgment for numerical error
INVALID_JUDGE = 1.0e30
INVALID_VALUE = 1.0e31
ZERO_JUDGE    = 1.0e-6
CUTOFF_UNDER_EXP = -50.0

  # judgement for warning
  # Output warning if bond angle is larger than WARN_ANGLE.
WARN_ANGLE = F_PI - 0.17e0  # Minimum angle [rad]
  # Output warning if bond length is longer than WARN_BOND.
WARN_BOND        = 5.0e0 # Maximum bond length [angst.]
WARN_BOND_RNA_SP = 6.0e0 # Maximum length for RNA S-P bond. [angst.]
WARN_BOND_RNA_SB = 8.0e0 # Maximum length for RNA S-B bond. [angst.]
WARN_BOND_RNA_PS = 6.0e0 # Maximum length for RNA P-S bond. [angst.]
WARN_RNA_O3_P = 1.8 # Maximum atomic distance b/w O3 and P in RNA [angst.]
WARN_RNA_P_O5 = 1.8 # Maximum atomic distance b/w P and O5 in RNA [angst.]

  # Mass
MASS_P = 30.973761e0  # Phosphorus
MASS_O = 15.9994e0    # Oxygen
MASS_C = 12.0107e0    # Carbon
MASS_N = 14.0067e0    # Nitrogen
MASS_BR= 79.904e0     # Boron
MASS_F = 18.9984032e0 # Fluorine
MASS_S = 32.065e0     # Sulfur
MASS_PO4 = MASS_P + MASS_O * 4.0e0 # Phosphoric acid

 # Nose-Hoover parameter
MXCS = 5               # # of thermal particle
#  real(PREC), parameter :: CSMASS = 100.0e0  # Mass of thermal particle
  
# O4'+C4'+C3'+C2'+C1'
MASS_RING     = MASS_C * 4 + MASS_O * 1
MASS_RIBOSE   = MASS_C * 5 + MASS_O * 2
MASS_ADENINE  = MASS_N * 5 + MASS_C * 5

MASS_GUANINE  = MASS_N * 5 + MASS_C * 5 + MASS_O * 1
MASS_URACIL   = MASS_N * 2 + MASS_C * 4 + MASS_O * 2
MASS_CYTOSINE = MASS_N * 3 + MASS_C * 4 + MASS_O * 1

MASS_RIN_SUGAR= MASS_P * 1 + MASS_C * 10 + MASS_O * 9
MASS_SUGAR    = MASS_C * 5 + MASS_O * 3
MASS_PHOS     = MASS_P * 1 + MASS_O * 2
