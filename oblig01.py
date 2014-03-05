import numpy as np

# Physical constants:
c     = 2.998e8       # [m/s]
sigma = 5.67e-8       # [W/m**2/K**4]
a     = 4 * sigma / c # [J/m**3/K**4]
u     = 1.660538921e-27 # [kg]

# Initial physical parameters for star:
L0   = 3.846e26       # [w] # L_sun
R0   = 0.5 * 6.96e8   # [m] # 0.5 * R_sun
M0   = 0.7 * 1.989e30 # [kg] # 0.7 * M_sun
rho0 = 1e3            # [kg/m**3]
T0   = 1e5            # [K]
P0   = 1e11           # [Pa]
# TODO Do no set all of rho0, T0 and P0. Calculate one of them.


# Read opacity:
infile = open("opacity.txt", "r")
logR_list = infile.readline()
logR_list = logR_list.split()[1:]
infile.readline()
logT_list = []
kappa_table = []
for line in infile:
    line = line.split()
    logT_list.append(line.pop(0))
    kappa_table.append(line)

logT_list, logR_list, kappa_table = np.array(logT_list), np.array(logR_list), np.array(kappa_table)
logT_list   = logT_list.astype(float)
logR_list   = logR_list.astype(float)
kappa_table = kappa_table.astype(float)


# Ratios for star:
X     = 0.7
Y_3   = 1e-10
Y     = 0.29
Z     = 0.01
Z_7Li = 1e-5
Z_7Be = 1e-5

# Numerical parameters:
n = 10000

# Set arrays:
L      = np.zeros(n+1)
L[0]   = L0
R      = np.zeros(n+1)
R[0]   = R0
M      = np.zeros(n+1)
M[0]   = M0
rho    = np.zeros(n+1)
rho[0] = rho0
T      = np.zeros(n+1)
T[0]   = T0
P      = np.zeros(n+1)
P[0]   = P0


# Particle list:
class Particle:
    u = 1.660538921e-27 # [kg] # atomic mass unit

    def __init__(s, mass, ratio):
        if mass > 0.5:
            mass *= u # assume mass was given in [u] and convert to [kg]
    
        s.mass  = mass  # particle mass
        s.ratio = ratio # total mass ratio

H   = Particle(1.6738e-27, X)
He3 = Particle(5.0081e-27, Y_3)
He4 = Particle(6.6464e-27, Y - Y_3)
Li7 = Particle(7.01600455, Z_7Li)
Be7 = Particle(7.01692983, Z_7Be)


# Functions:
def lam(i, j, T):
    T = T / 1e9

    if i == H and j == H:
        return 4.01e-15 * T**(-2/3.) * np.exp(- 3.380 * T**(-1/3.)) \
               * ( 1 + 0.123 * T**(1/3.) + 1.09 * T**(2/3.) + 0.938 * T )
    
    #if i == He3 and j == H3:
     #   return 6.04e10 * T**(-2/3.) * exp(- 12.276 T**(-1/3.)) \
      #         * ( 1 + 0.034 * T**(1/3.) - 0.522 * T**(2/3.) - 0.124 * T \
       #           + 0.353 * T**(4/3.) + 0.213 * T**(-5/3.) )
    
    if (i == H3 and j == H4) or (i == H4 and j == H3):
        T_star = T / (1 + 4.95e-2 * T)
        return 5.61e6 * T_star**(5/6.) * T**(-3/2.) \
               * exp(- 12.826 * T_star**(-1/3.))


def r(i, j, rho):

    if i == j:
        delta = 1
    else:
        delta = 0

    n_i = i.ratio * rho / i.mass
    n_j = j.ratio * rho / j.mass
    return n_i * n_j / (rho * (1 + delta))


def kappa(T, rho):

    logT = np.log10(T)
    rho = rho * 0.001 # convert to [g/cm**3]
    logR = np.log10(rho / T * 1e-6)

    i = 0 # this will be the vertical index in kappa_table
    for logT_l in logT_list:
        if logT_l > logT: # compare sizes to find the closes
            if i > 0:
                if logT - logT_list[i-1] < logT_l - logT: # check if the previous value was closer
                    i -= 1
            break # stop checking
        i += 1
        
    j = 0 # this will be the horizontal index in kappa_table
    for logR_l in logR_list:
        if logR_l > logR:
            if j > 0:
                if logR - logR_list[j-1] < logR_l - logR:
                    j -= 1
            break
        j += 1
    
    kappa = 10**kappa_table[i,j] # find value in table and convert back from log
    return kappa * 1000 # convert to [kg/m**3]

    
"""
# Integration loop:
for i in range(n):

    eps = 

    r[i+1] = 1 / (4 * np.pi * r[i]**2 * rho[i]) * dm
    P[i+1] = - G * m[i] / (4 * np.pi * r[i]) * dm
    L[i+1] = eps * dm
    P_rad[i+1] = a / 3. * T[i]**4
    P_gas[i+1] = P - P_rad[i]
    T[i+1] = - 3 * kappa(T[i], rho[i]) * L[i] / (256 * np.pi*np.pi * sigma * r[i]**2 * T[i]**3)
"""
