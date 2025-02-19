System  = NiO

KPAR = 5

# better convergence for small kpt grids
ISMEAR = 2
SIGMA = 0.05

# converge wave functions
EDIFF = 1.E-7
NELMIN = 35

# the energy window to optimize projector channels (absolute)
EMIN = -3
EMAX = 10

LMAXMIX = 6

# project to Ni d and O p states
LORBIT = 14
LOCPROJ = "1 : d : Pr
           2 : p : Pr"

# CSC flags
ICHARG = 5
NELM = 1000
NELMIN = 1000
NELMDL = -1
IMIX=1
BMIX=0.5
AMIX=0.02
LSYNCH5=True
