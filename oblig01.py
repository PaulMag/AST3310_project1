import numpy as np

# Physical constants:
c     = 2.998e8       # [m/s]
sigma = 5.67e-8       # [W/m**2/K**4]
a     = 4 * sigma / c # [J/m**3/K**4]

# Initial physical parameters for star:
L0   = 3.846e26       # [w] # L_sun
R0   = 0.5 * 6.96e8   # [m] # 0.5 * R_sun
M0   = 0.7 * 1.989e30 # [kg] # 0.7 * M_sun
rho0 = 1e3            # [kg/m**3]
T0   = 1e5            # [K]
P0   = 1e11           # [Pa]
# TODO Do no set all of rho0, T0 and P0. Calculate one of them.
# TODO Import opacities kappa.

# Ratios for star:
X     = 0.7
Y_3   = 1e-10
Y     = 0.29
Z     = 0.01
Z_7Li = 1e-5
Z_7Be = 1e-5

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


# Numerical parameters:
n = 10000

# Functions:
T9 = lambda T = T / 1e9

# Integration loop:
for i in range(n):

    r[i+1] = 1 / (4 * np.pi * r[i]**2 * rho[i]) * dm
    P[i+1] = - G * m[i] / (4 * np.pi * r[i]) * dm
    L[i+1] = eps * dm
    P_rad[i+1] = a / 3. * T[i]**4
    P_gas[i+1] = P - P_rad
    T[i+1] = - 3 * kappa[i] * L[i] / (256 * np.pi*np.pi * sigma * r[i]**2 * T[i]**3)
