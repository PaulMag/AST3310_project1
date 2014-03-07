import numpy as np

# Physical constants:
c     = 2.998e8         # [m/s]
sigma = 5.67e-8         # [W/m**2/K**4]
a     = 4 * sigma / c   # [J/m**3/K**4]
u     = 1.660538921e-27 # [kg]
k     = 1.3806488e-23   # [J/K]
G     = 6.67384e-11     # [m**3/kg/s**2]
MeV   = 1.602176565e-13 # [J]
avogadro_inverse = 1 / 6.0221413e23

# Initial physical parameters for star:
L0   = 3.846e26       # [w] # L_sun
R0   = 0.5 * 6.96e8   # [m] # 0.5 * R_sun
M0   = 0.7 * 1.989e30 # [kg] # 0.7 * M_sun
rho0 = 1e3            # [kg/m**3]
T0   = 1e5            # [K]
P0   = 1e11           # [Pa]

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
infile.close()

logT_list, logR_list, kappa_table = np.array(logT_list), np.array(logR_list), np.array(kappa_table)
logT_list   = logT_list.astype(float)
logR_list   = logR_list.astype(float)
kappa_table = kappa_table.astype(float)
kappa_x = len(logR_list)
kappa_y = len(logT_list)


# Ratios for star:
X     = 0.7
Y_3   = 1e-10
Y     = 0.29
Z     = 0.01
Z_7Li = 1e-5
Z_7Be = 1e-5

mu = 1. / (2*X + 3*Y/4. + Z/2.)


# Q's:
Q_p_p     =  1.177 * MeV
Q_d_p     =  5.494 * MeV
Q_He3_He3 = 12.860 * MeV

Q_He3_He4 =  1.586 * MeV
Q_Be7_e   =  0.049 * MeV
Q_Li7_p   = 17.346 * MeV

Q_Be7_p   =  0.137 * MeV
Q_B8      =  8.367 * MeV
Q_Be8     =  2.995 * MeV

Q_123 = Q_p_p + Q_d_p # the numbers after Q_ tells which PP chains used this Q
Q_1   = Q_He3_He3
Q_23  = Q_He3_He4
#Q_2   = Q_Be7_e + Q_Li7_p      # can I add these? # ended up not
Q_2a  = Q_Be7_e
Q_2b  = Q_Li7_p
Q_3   = Q_Be7_p + Q_B8 + Q_Be8 # can I add these? # Boris said yes, since missing

# Numerical parameters:
n = int(1e12)
dm = - M0 / float(n)
outfile = open("luminosity.dat", "w")


# Set arrays:
array_size = 1 # + 1

L      = np.zeros(array_size+1)
L[0]   = L0

R      = np.zeros(array_size+1)
R[0]   = R0

M      = np.zeros(array_size+1)
M[0]   = M0

rho    = np.zeros(array_size+1)
rho[0] = rho0

T      = np.zeros(array_size+1)
T[0]   = T0

P      = np.zeros(array_size+1)
P[0]   = P0

P_rad    = np.zeros(array_size+1)
P_rad[0] = a / 3. * T[0]**4

P_gas    = np.zeros(array_size+1)
P_gas[0] = P[0] - P_rad[0]

# Do no set all of rho0, T0 and P0. Calculate one of them:
rho[0] = P_gas[0] * mu * u / (k * T[0])


# Particle list:
class Particle:

    def __init__(s, mass, ratio):
        # particle mass, particle ratio, relative particle density
        # n = rho_rel * rho
        if mass > 0.5:
            mass *= u # assume mass was given in [u] and convert to [kg]
    
        s.mass  = mass  # particle mass
        s.ratio = ratio # total mass ratio
        s.rho_rel = ratio / mass

H   = Particle(1.6738e-27, X)
He3 = Particle(5.0081e-27, Y_3)
He4 = Particle(6.6464e-27, Y - Y_3)
Li7 = Particle(7.01600455, Z_7Li)
Be7 = Particle(7.01692983, Z_7Be)
# Make the electron and set relative particle density to n_e = n_H + n_He:
e_  = Particle(9.10938291e-31, 0)
e_.rho_rel = H.rho_rel + He3.rho_rel + He4.rho_rel

# Functions:
def lam(i, j, T):
    T = T / 1e9

    # PP I & II & III:
    if i == H and j == H:
        s = 4.01e-15 * T**(-2/3.) * np.exp(- 3.380 * T**(-1/3.)) \
            * ( 1 + 0.123 * T**(1/3.) + 1.09 * T**(2/3.) + 0.938 * T )
    
    # PP I:
    if i == He3 and j == He3:
        s = 6.04e10 * T**(-2/3.) * np.exp(- 12.276 * T**(-1/3.)) \
            * ( 1 + 0.034 * T**(1/3.) - 0.522 * T**(2/3.) - 0.124 * T \
               + 0.353 * T**(4/3.) + 0.213 * T**(-5/3.) )
    
    # PP II & III:
    if (i == He3 and j == He4) or (i == He4 and j == He3):
        T_star = T / (1 + 4.95e-2 * T)
        s = 5.61e6 * T_star**(5/6.) * T**(-3/2.) \
            * np.exp(- 12.826 * T_star**(-1/3.))

    # PP II:
    if (i == Be7 and j == e_) or (i == e_ and j == Be7):
        if T * 1e3 < 1: # T < 1e6
            s = 1.57e-7 / (e_.rho_rel * rho[0])
        else:
            s = 1.34e-10 * T**(-1/2.) * (1 - 0.537 * T**(1/3.) + 3.86 * T**(2/3.)
                + 0.0027 * T**(-1) * np.exp(2.515e-3 * T**(-1)))
    
    if (i == Li7 and j == H) or (i == H and j == Li7):
        T_star = T / (1 + 0.759 * T)
        s = 1.096e9 * T**(-2/3.) * np.exp(- 8.472 * T**(-1/3.)) \
            - 4.830e8 * T_star**(5/6.) * T**(-2/3.) \
            * np.exp(- 8.472 * T_star**(-1/3.)) \
            + 1.06e10 * T**(-3/2.) * np.exp(- 30.442 * T**(-1))

    # PP III:
    if (i == Be7 and j == H) or (i == H and j == Be7):
        s = 3.11e5 * T**(-2/3.) * np.exp(- 10.262 * T**(-1/3.)) \
            + 2.53e3 * T**(-3/2.) * np.exp(- 7.306 * T**(-1))

    #print "LAMBDA = ", s * avogadro_inverse * 1e-6 # for debugging
    return s * 1e-6 # convert to [m**3/s]


def rate(i, j, rho, T):

    if i == j:
        delta_1 = 2 # this is delta + 1
    else:
        delta_1 = 1

    # n_j * n_i / (rho + (1 + delta)) * lambda_ij:
    return i.rho_rel * j.rho_rel * rho / delta_1 * lam(i, j, T)


def kappa(T, rho):

    logT = np.log10(T)
    rho = rho * 0.001 # convert to [g/cm**3]
    logR = np.log10(rho / T * 1e-6)

    i = 0 # this will be the vertical index in kappa_table
    for logT_l in logT_list:
        if logT_l > logT: # compare sizes to find the closest
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
    
    if i >= kappa_y:
        i = kappa_y - 1
        print "Warning: T may have blown the kappa table."
    
    if j >= kappa_x:
        j = kappa_x - 1
        print "Warning: rho may have blown the kappa table."
    
    kappa = 10**kappa_table[i,j] # find value in table and convert back from log
    return kappa * 1000 # convert to [kg/m**3]


outfile.write(str(L[0]))
print rho[0], R[0], P[0], L[0], P_rad[0], P_gas[0], T[0]
print "dm = ", dm
# Integration loop:
for i in range(n / 1000000000):

    eps =   rate(H  , H,   rho[0], T[0]) * Q_123 \
          + rate(He3, He3, rho[0], T[0]) * Q_1   \
          + rate(He3, He4, rho[0], T[0]) * Q_23  \
          + rate(Be7, e_,  rho[0], T[0]) * Q_2a  \
          + rate(Li7, H,   rho[0], T[0]) * Q_2b  \
          + rate(Be7, H,   rho[0], T[0]) * Q_3
    eps *=  avogadro_inverse
    print "eps = ", eps # for debugging
    
    rho[1] = P_gas[0] * mu * u / (k * T[0])
    P_rad[1] = a / 3. * T[0]**4
    P_gas[1] = P[0] - P_rad[0]

    R[1] = R[0] + 1. / (4 * np.pi * R[0]**2 * rho[0]) * dm
    P[1] = P[0] - G * M[0] / (4 * np.pi * R[0]) * dm
    L[1] = L[0] + eps * dm
    T[1] = T[0] - 3 * kappa(T[0], rho[0]) * L[0] \
                  / (256 * np.pi*np.pi * sigma * R[0]**4 * T[0]**3) * dm

    # Reset stuff:
    rho[0] = rho[1]
    R[0]   = R[1]
    P[0]   = P[1]
    L[0]   = L[1]
    P_rad[0] = P_rad[1]
    P_gas[0] = P_gas[1]
    T[0]   = T[1]

    # Check if anything dropped below zero:
    if rho[0] <= 0 or R[0] <= 0 or P[0] <= 0 or L[0] <= 0 or P_rad[0] <= 0 \
       or P_gas[0] <= 0 or T[0] <= 0:
        print "Something dropped below 0. Stop simulation."
        print "Progress = %d / %d" % (i, n)
        print rho[0], R[0], P[0], L[0], P_rad[0], P_gas[0], T[0]
        break

    # Write to file:
    outfile.write(str(L[0]) + "\n")

